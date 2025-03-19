#include <inttypes.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <assert.h>
#include <getopt.h>
#include <err.h>
#include <assert.h>

typedef uint64_t u64;       /* portable 64-bit integer */
typedef uint32_t u32;       /* portable 32-bit integer */
struct __attribute__ ((packed)) entry { u32 k; u64 v; };  /* hash table entry */
struct __attribute__ ((packed)) write_pair { u64 f_x; u64 x; };

/***************************** global variables ******************************/

u64 n = 0;         /* block size (in bits) */
u64 mask;          /* this is 2**n - 1 */

int buffer_size = 1024;

u64 global_dict_size;     /* number of slots in the global hash table */
u64 local_dict_size;



u64 min_local_h;
u64 max_local_h;
struct entry *A;   /* the hash table */
struct write_pair *write_buffer;
int *write_buffer_size;
u64 *read_buffer;
int *read_buffer_size;

void *send_write_buffer;
void *receive_write_buffer;
void *send_read_buffer;
void *receive_read_buffer;
void *send_results_buffer;
void *receive_results_buffer;

MPI_Request *send_requests;
MPI_Request *receive_requests;

int *finished;

/* (P, C) : two plaintext-ciphertext pairs */
u32 P[2][2] = {{0, 0}, {0xffffffff, 0xffffffff}};
u32 C[2][2];

int world_rank, world_size;

/************************ tools and utility functions *************************/

double wtime()
{
	struct timeval ts;
	gettimeofday(&ts, NULL);
	return (double)ts.tv_sec + ts.tv_usec / 1E6;
}

// murmur64 hash functions, tailorized for 64-bit ints / Cf. Daniel Lemire
u64 murmur64(u64 x)
{
    x ^= x >> 33;
    x *= 0xff51afd7ed558ccdull;
    x ^= x >> 33;
    x *= 0xc4ceb9fe1a85ec53ull;
    x ^= x >> 33;
    return x;
}

/* represent n in 4 bytes */
void human_format(u64 n, char *target)
{
    if (n < 1000) {
        sprintf(target, "%" PRId64, n);
        return;
    }
    if (n < 1000000) {
        sprintf(target, "%.1fK", n / 1e3);
        return;
    }
    if (n < 1000000000) {
        sprintf(target, "%.1fM", n / 1e6);
        return;
    }
    if (n < 1000000000000ll) {
        sprintf(target, "%.1fG", n / 1e9);
        return;
    }
    if (n < 1000000000000000ll) {
        sprintf(target, "%.1fT", n / 1e12);
        return;
    }
}

/******************************** SPECK block cipher **************************/

#define ROTL32(x,r) (((x)<<(r)) | (x>>(32-(r))))
#define ROTR32(x,r) (((x)>>(r)) | ((x)<<(32-(r))))

#define ER32(x,y,k) (x=ROTR32(x,8), x+=y, x^=k, y=ROTL32(y,3), y^=x)
#define DR32(x,y,k) (y^=x, y=ROTR32(y,3), x^=k, x-=y, x=ROTL32(x,8))

void Speck64128KeySchedule(const u32 K[],u32 rk[])
{
    u32 i,D=K[3],C=K[2],B=K[1],A=K[0];
    for(i=0;i<27;){
        rk[i]=A; ER32(B,A,i++);
        rk[i]=A; ER32(C,A,i++);
        rk[i]=A; ER32(D,A,i++);
    }
}

void Speck64128Encrypt(const u32 Pt[], u32 Ct[], const u32 rk[])
{
    u32 i;
    Ct[0]=Pt[0]; Ct[1]=Pt[1];
    for(i=0;i<27;)
        ER32(Ct[1],Ct[0],rk[i++]);
}

void Speck64128Decrypt(u32 Pt[], const u32 Ct[], u32 const rk[])
{
    int i;
    Pt[0]=Ct[0]; Pt[1]=Ct[1];
    for(i=26;i>=0;)
        DR32(Pt[1],Pt[0],rk[i--]);
}

/******************************** dictionary ********************************/

/*
 * "classic" hash table for 64-bit key-value pairs, with linear probing.
 * It operates under the assumption that the keys are somewhat random 64-bit integers.
 * The keys are only stored modulo 2**32 - 5 (a prime number), and this can lead
 * to some false positives.
 */
static const u32 EMPTY = 0xffffffff;
static const u64 PRIME = 0xfffffffb;

void dict_setup(u64 size);
void initialise_buffers(int size);
int store_write_in_buffer(u64 f_x, u64 x);
int store_read_in_buffer(u64 g_y);
void dict_insert_sender(int procNo);
void dict_insert_receiver(int procNo);
int dict_probe_sender(u64 g_y[], int size, int maxval, int nval[], u64 values[]);
void dict_probe_receiver(u64 g_y[], int size, int maxval, int nval[], u64 values[]);


/* allocate a hash table with `size` slots (12*size bytes) */
void dict_setup(u64 size)
{
	global_dict_size = size;
    local_dict_size = size / world_size;
    min_local_h = local_dict_size * world_rank;
    max_local_h = local_dict_size * (world_rank + 1);

	char hdsize[8];
	human_format(local_dict_size * sizeof(*A), hdsize);
	printf("Rank %d Dictionary size: %sB\n", world_rank, hdsize);

	A = (struct entry *) malloc(sizeof(*A) * local_dict_size);
	if (A == NULL)
		err(1, "impossible to allocate the dictionnary");
    // TODO: parallelise with OpenMP
	for (u64 i = 0; i < local_dict_size; i++)
		A[i].k = EMPTY;
}

void initialise_buffers(int size) {
    write_buffer = (struct write_pair *) malloc(sizeof(*write_buffer) * size * world_size);
    write_buffer_size = (int *) malloc(sizeof(int) * world_size);
    read_buffer = (u64 *) malloc(sizeof(*read_buffer) * size * world_size);
    read_buffer_size = (int *) malloc(sizeof(int) * world_size);

    send_write_buffer = (void *) malloc((sizeof(int) + sizeof(struct write_pair) * size) * world_size);
    receive_write_buffer = (void *) malloc((sizeof(int) + sizeof(struct write_pair) * size) * world_size);

    send_requests = (MPI_Request *) malloc(sizeof(MPI_Request) * world_size);
    receive_requests = (MPI_Request *) malloc(sizeof(MPI_Request) * world_size);

    finished = malloc(sizeof(int) * world_size);

    for (int i = 0; i < size; ++i) {
        write_buffer_size[i] = read_buffer_size[i] = 0;
    }

    for (int i = 0; i < world_size; ++i) {
        finished[i] = 0;
    }
}

/* Stores a f(x) x pair in the buffer to be sent */
int store_write_in_buffer(u64 f_x, u64 x) {
    u64 h = murmur64(f_x) % global_dict_size;
    int proc = h / local_dict_size;

    write_buffer[proc * buffer_size + write_buffer_size[proc]].f_x = f_x;
    write_buffer[proc * buffer_size + write_buffer_size[proc]].x = x;
    write_buffer_size[proc]++;
    if (write_buffer_size[proc] == buffer_size) {
        dict_insert_sender(proc);
        write_buffer_size[proc] = 0;
    }
}

/* Stores a g(y) entry in the buffer to be sent */
int store_read_in_buffer(u64 g_y) {
    u64 h = murmur64(g_y) % global_dict_size;
    int proc = h / local_dict_size;

    read_buffer[proc * buffer_size + read_buffer_size[proc]] = g_y;
    read_buffer_size[proc]++;
    if (read_buffer_size[proc] == buffer_size) {
        int *nval = malloc(sizeof(int) * buffer_size);
        u64 *values = malloc(sizeof(u64) * buffer_size * 256);
        dict_probe_sender(read_buffer,buffer_size,256,nval,values);
    }
}

void dict_insert_sender(int procNo) {
    // MPI_Send
    struct write_pair *local_write_buffer = write_buffer + procNo * buffer_size;
    int num_entries = write_buffer_size[procNo];

    if (num_entries == 0)
        return;

    int swb_cap = sizeof(int) + buffer_size * sizeof(struct write_pair);
    void *ptr = send_write_buffer + (swb_cap * procNo);
    int *int_ptr = (int *) ptr;
    *int_ptr = num_entries;
    struct write_pair *wp_ptr = (struct write_pair *)(int_ptr+1);
    // TODO: OpenMPI
    for (int i = 0; i < num_entries; ++i) {
        wp_ptr[i] = local_write_buffer[i];
    }
    write_buffer_size[procNo] = 0;
    MPI_Isend(ptr, swb_cap, MPI_BYTE, procNo, 0, MPI_COMM_WORLD, send_requests + procNo);
}

void dict_insert_receiver(int procNo) {
    int swb_cap = sizeof(int) + buffer_size * sizeof(struct write_pair);
    void *ptr = receive_write_buffer + swb_cap * procNo;
    int *int_ptr = (int *) ptr;
    int num_entries = *int_ptr;
    struct write_pair *wp_ptr = (struct write_pair *) (int_ptr + 1);

    if (num_entries == 0) {
        finished[procNo] = 1;
        return;
    }

    // TODO: parallelise with OpenMP
    for (int i = 0; i < num_entries; ++i) {
        u64 h = murmur64(wp_ptr[i].f_x) % global_dict_size;
        u32 h_p = wp_ptr[i].f_x % PRIME;
        u64 x = wp_ptr[i].x;
        assert(h >= min_local_h && h <= max_local_h);
        u64 local_h = h % local_dict_size;
        assert(local_h >= min_local_h && local_h <= max_local_h);
        u64 insertion_point = local_h;
        while (A[insertion_point].k != EMPTY) {
            assert((++insertion_point) % local_dict_size != local_h);
        }
        A[insertion_point].k = h_p;
        A[insertion_point].v = x;
    }

    MPI_Irecv(ptr, swb_cap, MPI_BYTE, procNo, 0, MPI_COMM_WORLD, receive_requests + procNo);
}

/* Insert the binding key |----> value in the dictionnary */
/*void dict_insert(u64 key, u64 value)
{
    u64 h = murmur64(key) % dict_size;
    for (;;) {
        if (A[h].k == EMPTY)
            break;
        h += 1;
        if (h == dict_size)
            h = 0;
    }
    assert(A[h].k == EMPTY);
    A[h].k = key % PRIME;
    A[h].v = value;
}*/

int dict_probe_sender(u64 g_y[], int size, int maxval, int nval[], u64 values[]) {
    // MPI_Send
}

/* Query the dictionnary with this `key`.  Write values (potentially)
 *  matching the key in `values` and return their number. The `values`
 *  array must be preallocated of size (at least) `maxval`.
 *  The function returns -1 if there are more than `maxval` results.
 */
void dict_probe_receiver(u64 g_y[], int size, int maxval, int nval[], u64 values[])
{
    // TODO: OpenMP
    for (int i = 0; i < size; ++i) {
        u64 h = murmur64(g_y[i]) % global_dict_size;
        u32 h_p = g_y[i] % PRIME;
        u64 local_h = h % local_dict_size;
        assert(h >= min_local_h && h <= max_local_h);

        nval[i] = 0;
        for (;;h = (h+1)%local_dict_size) {
            if (A[h].k == EMPTY)
                break;
            if (A[h].k == h_p) {
                if (nval[i] == maxval) {
                    nval[i] = -1;
                    break;
                }
                values[i*maxval + nval[i]] = A[h].v;
                nval[i]++;
            }
        }
    }
}

/***************************** MITM problem ***********************************/

/* f : {0, 1}^n --> {0, 1}^n.  Speck64-128 encryption of P[0], using k */
u64 f(u64 k)
{
    assert((k & mask) == k);
    u32 K[4] = {k & 0xffffffff, k >> 32, 0, 0};
    u32 rk[27];
    Speck64128KeySchedule(K, rk);
    u32 Ct[2];
    Speck64128Encrypt(P[0], Ct, rk);
    return ((u64) Ct[0] ^ ((u64) Ct[1] << 32)) & mask;
}

/* g : {0, 1}^n --> {0, 1}^n.  speck64-128 decryption of C[0], using k */
u64 g(u64 k)
{
    assert((k & mask) == k);
    u32 K[4] = {k & 0xffffffff, k >> 32, 0, 0};
    u32 rk[27];
    Speck64128KeySchedule(K, rk);
    u32 Pt[2];
    Speck64128Decrypt(Pt, C[0], rk);
    return ((u64) Pt[0] ^ ((u64) Pt[1] << 32)) & mask;
}

bool is_good_pair(u64 k1, u64 k2)
{
    u32 Ka[4] = {k1 & 0xffffffff, k1 >> 32, 0, 0};
    u32 Kb[4] = {k2 & 0xffffffff, k2 >> 32, 0, 0};
    u32 rka[27];
    u32 rkb[27];
    Speck64128KeySchedule(Ka, rka);
    Speck64128KeySchedule(Kb, rkb);
    u32 mid[2];
    u32 Ct[2];
    Speck64128Encrypt(P[1], mid, rka);
    Speck64128Encrypt(mid, Ct, rkb);
    return (Ct[0] == C[1][0]) && (Ct[1] == C[1][1]);
}

/******************************************************************************/

/* search the "golden collision" */
int golden_claw_search(int maxres, u64 k1[], u64 k2[])
{
    double start = wtime();
    u64 N = 1ull << n;
    // setup receive request listeners
    int sbw_cap = (sizeof(int) + sizeof(struct write_pair) * buffer_size);
    for (int i = 0; i < world_size; ++i) {
        if (i == world_rank) continue;
        MPI_Irecv(receive_write_buffer + sbw_cap * i, sbw_cap, MPI_BYTE, i, 0, MPI_COMM_WORLD, receive_requests + i);
    }
    // TODO: OpenMP
    for (u64 x = (world_rank * N) / world_size; x < ((world_rank+1) * N) / world_size; ++x) {
        u64 z = f(x);
        store_write_in_buffer(z, x);
        if (x % 100 == 0) { // check receive requests
            int req_recvd;
            for (int i = 0; i < world_size; ++i) {
                if (i == world_rank || finished[i]) continue;
                MPI_Test(receive_requests + i, &req_recvd, MPI_STATUS_IGNORE);
                if (req_recvd) {
                    dict_insert_receiver(i);
                }
            }
        }
    }

    for (int i = 0; i < world_size; ++i) {
        if (i == world_rank) continue;
        if (write_buffer_size[i] != 0) { // flush buffers
            dict_insert_sender(i);
        }
        dict_insert_sender(i); // signal to terminate listeners
    }

    int fw = 0;
    while (!fw) {
        fw = 1;
        for (int i = 0; i < world_size; ++i) {
            if (i == world_rank || finished[i]) continue;
            int req_recvd;
            MPI_Test(receive_requests + i, &req_recvd, MPI_STATUS_IGNORE);
            if (req_recvd) {
                dict_insert_receiver(i);
                if (!finished[i])
                    fw = 0;
            } else {
                fw = 0;
            }
        }
    }

    double mid = wtime();
    printf("Fill: %.1fs\n", mid - start);



    int nres = 0;

    for (int i = 0; i < world_size; ++i) {
        if (i == world_rank) continue;
        MPI_Irecv(receive_read_buffer +
    }

    for (int i = 0; i < world_size; ++i) {
        finished[i] = 0;
    }

    for (u64 y = (world_rank * N) / world_size; y < ((world_rank+1) * N) / world_size; ++y) {
        u64 g_y = g(y);
        store_read_in_buffer(g_y,y);
        if (y % 100 == 0) {
            int req_recvd;
            for (int i = 0; i < world_size; ++i) {
                if (i == world_rank || finished[i]) continue;
                MPI_Test(receive_requests + i, &req_recvd, MPI_STATUS_IGNORE);
            }
        }
    }

    /*int nres = 0;
    u64 ncandidates = 0;
    u64 x[256];
    for (u64 z = 0; z < N; z++) {
        u64 y = g(z);

        store_read_in_buffer(y);

        assert(nx >= 0);
        ncandidates += nx;
        for (int i = 0; i < nx; i++)
            if (is_good_pair(x[i], z)) {
            	if (nres == maxres)
            		return -1;
            	k1[nres] = x[i];
            	k2[nres] = z;
            	printf("SOLUTION FOUND!\n");
            	nres += 1;
            }
    }
    printf("Probe: %.1fs. %" PRId64 " candidate pairs tested\n", wtime() - mid, ncandidates);
    return nres;*/
    return 0; // TODO: change
}

/************************** command-line options ****************************/

void usage(char **argv)
{
        printf("%s [OPTIONS]\n\n", argv[0]);
        printf("Options:\n");
        printf("--n N                       block size [default 24]\n");
        printf("--C0 N                      1st ciphertext (in hex)\n");
        printf("--C1 N                      2nd ciphertext (in hex)\n");
        printf("\n");
        printf("All arguments are required\n");
        exit(0);
}

void process_command_line_options(int argc, char ** argv)
{
        struct option longopts[4] = {
                {"n", required_argument, NULL, 'n'},
                {"C0", required_argument, NULL, '0'},
                {"C1", required_argument, NULL, '1'},
                {NULL, 0, NULL, 0}
        };
        char ch;
        int set = 0;
        while ((ch = getopt_long(argc, argv, "", longopts, NULL)) != -1) {
                switch (ch) {
                case 'n':
                        n = atoi(optarg);
                        mask = (1 << n) - 1;
                        break;
                case '0':
                        set |= 1;
                        u64 c0 = strtoull(optarg, NULL, 16);
                        C[0][0] = c0 & 0xffffffff;
                        C[0][1] = c0 >> 32;
                        break;
                case '1':
                        set |= 2;
                        u64 c1 = strtoull(optarg, NULL, 16);
                        C[1][0] = c1 & 0xffffffff;
                        C[1][1] = c1 >> 32;
                        break;
                default:
                        errx(1, "Unknown option\n");
                }
        }
        if (n == 0 || set != 3) {
        	usage(argv);
        	exit(1);
        }
}

/******************************************************************************/

int main(int argc, char **argv)
{
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	process_command_line_options(argc, argv);
    printf("Running with n=%d, C0=(%08x, %08x) and C1=(%08x, %08x)\n",
        (int) n, C[0][0], C[0][1], C[1][0], C[1][1]);

	dict_setup(1.125 * (1ull << n));
    initialise_buffers(buffer_size);

	/* search */
	u64 k1[16], k2[16];
	int nkey = golden_claw_search(16, k1, k2);
	assert(nkey > 0);

	/* validation */
	for (int i = 0; i < nkey; i++) {
    	assert(f(k1[i]) == g(k2[i]));
    	assert(is_good_pair(k1[i], k2[i]));
	    printf("Solution found: (%" PRIx64 ", %" PRIx64 ") [checked OK]\n", k1[i], k2[i]);
	}

}