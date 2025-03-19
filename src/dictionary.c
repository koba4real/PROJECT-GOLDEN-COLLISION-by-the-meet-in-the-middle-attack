#include <inttypes.h>
#include <stdlib.h>
#include <err.h>
#include <string.h>
#include "dictionary.h"
#include "util.h"


/* tags:
 * 0 - write size send
 * 1 - write data send
 * 2 - read size send
 * 3 - read data send
 */

// global variables

struct entry *local_dictionary;
MPI_Comm world;
int world_size, world_rank, num_solutions;
u64 local_dict_size, global_dict_size, min_local_h, max_local_h;
struct write_request *write_send_buffer;
int *write_send_buffer_size;
struct write_request *write_recv_buffer;
int *write_recv_buffer_size;
struct read_request *read_send_buffer;
int *read_send_buffer_size;
struct read_request *read_recv_buffer;
int *read_recv_buffer_size;
MPI_Request *write_listeners;
MPI_Request *write_size_speakers;
MPI_Request *write_data_speakers;
MPI_Request *read_listeners;
MPI_Request *read_size_speakers;
MPI_Request *read_data_speakers;
int maxres;
int *open_active_write;
int *open_active_read;
u64 local_candidates;
struct write_request *write_recv_buffer;
struct read_request *read_recv_buffer;
struct result solutions[MAX_RESULTS];
/*omp_lock_t *write_send_locks;
omp_lock_t write_recv_lock;
omp_lock_t write_check_lock;
omp_lock_t *read_send_locks;
omp_lock_t read_recv_lock;
omp_lock_t read_check_lock;
omp_lock_t sol_lock;*/

void setup_dict(u64 dict_size, MPI_Comm comm_world, int max_results) {
    global_dict_size = dict_size;
    world = comm_world;
    MPI_Comm_size(world, &world_size);
    MPI_Comm_rank(world, &world_rank);

    num_solutions = 0;
    local_candidates = 0;
    maxres = max_results;

    // calculate local dict size
    local_dict_size = (global_dict_size + world_size - 1) / world_size;
    min_local_h = local_dict_size * world_rank;
    max_local_h = local_dict_size * (world_rank + 1);
    if (global_dict_size < max_local_h)
        max_local_h = global_dict_size;

    local_dictionary = malloc(sizeof(struct entry) * local_dict_size);
    if (local_dictionary == NULL) {
        err(1, "impossible to allocate the dictionary");
    }

    #pragma omp parallel for
    for (u64 i = 0; i < local_dict_size; i++) {
        local_dictionary[i].k = EMPTY;
    }

    // setup buffers
    write_send_buffer = malloc(sizeof(struct write_request) * BUFFER_SIZE * world_size);
    write_send_buffer_size = malloc(sizeof(u32) * world_size);
    write_recv_buffer = malloc(sizeof(struct write_request) * BUFFER_SIZE * world_size);
    write_recv_buffer_size = malloc(sizeof(u32) * world_size);

    read_send_buffer = malloc(sizeof(struct read_request) * BUFFER_SIZE * world_size);
    read_send_buffer_size = malloc(sizeof(u32) * world_size);
    read_recv_buffer = malloc(sizeof(struct read_request) * BUFFER_SIZE * world_size);
    read_recv_buffer_size = malloc(sizeof(u32) * world_size);

    for (int i = 0; i < world_size; i++) {
        write_send_buffer_size[i] = 0;
        write_recv_buffer_size[i] = 0;
        read_send_buffer_size[i] = 0;
        read_recv_buffer_size[i] = 0;
    }

    // setup listeners
    write_listeners = malloc(sizeof(MPI_Request) * world_size);
    read_listeners = malloc(sizeof(MPI_Request) * world_size);
    write_size_speakers = malloc(sizeof(MPI_Request) * world_size);
    write_data_speakers = malloc(sizeof(MPI_Request) * world_size);
    read_size_speakers = malloc(sizeof(MPI_Request) * world_size);
    read_data_speakers = malloc(sizeof(MPI_Request) * world_size);

    open_active_write = malloc(sizeof(int) * world_size);
    open_active_read = malloc(sizeof(int) * world_size);

    /*write_send_locks = malloc(sizeof(omp_lock_t) * world_size);
    read_send_locks = malloc(sizeof(omp_lock_t) * world_size);
    omp_init_lock(&write_recv_lock);
    omp_init_lock(&write_check_lock);
    omp_init_lock(&read_recv_lock);
    omp_init_lock(&read_check_lock);*/

    for (int i = 0; i < world_size; i++) {
        if (i != world_rank) {
            MPI_Irecv(write_recv_buffer_size + i, 1, MPI_INT, i, 0, world, write_listeners + i);
            MPI_Irecv(read_recv_buffer_size + i, 1, MPI_INT, i, 2, world, read_listeners + i);
        }
        write_size_speakers[i] = MPI_REQUEST_NULL;
        write_data_speakers[i] = MPI_REQUEST_NULL;
        read_size_speakers[i] = MPI_REQUEST_NULL;
        read_data_speakers[i] = MPI_REQUEST_NULL;

        open_active_write[i] = 1 * (i != world_rank) - 1 * (i == world_rank);
        open_active_read[i] = 1 * (i != world_rank) - 1 * (i == world_rank);

        //omp_init_lock(&write_send_locks[i]);
        //omp_init_lock(&read_send_locks[i]);
    }
}

int write_to_dict(u64 f_x, u64 x) {
    u64 h_f_x = murmur64(f_x);
    u64 indx = h_f_x % global_dict_size;
    int proc_num = indx / local_dict_size;

    u32 p_f_x = f_x % PRIME;
    if (proc_num == world_rank) { // store locally
        u64 local_index = indx - min_local_h;
        //omp_set_lock(&write_recv_lock);
        while (local_dictionary[local_index].k != EMPTY) {
            local_index = (local_index + 1) % local_dict_size;
            if (local_index == indx - min_local_h) {
                fprintf(stderr, "Error: no space in local dictionary");
                exit(1);
            }
        }
        local_dictionary[local_index].k = p_f_x;
        local_dictionary[local_index].v = x;
        //omp_unset_lock(&write_recv_lock);
    } else {
        //omp_set_lock(&write_send_locks[proc_num]);
        if (write_send_buffer_size[proc_num] == BUFFER_SIZE) {
            int size_sent, data_sent;
            MPI_Test(write_size_speakers + proc_num, &size_sent, MPI_STATUS_IGNORE);
            MPI_Test(write_data_speakers + proc_num, &data_sent, MPI_STATUS_IGNORE);
            if (!size_sent || !data_sent) {
                //omp_unset_lock(&write_send_locks[proc_num]);
                return 0;
            }
            write_send_buffer_size[proc_num] = 0;
        }
        struct write_request *addr = write_send_buffer + (proc_num * BUFFER_SIZE + write_send_buffer_size[proc_num]);
        addr->f_x = f_x;
        addr->x = x;
        write_send_buffer_size[proc_num]++;

        if (write_send_buffer_size[proc_num] == BUFFER_SIZE) {
            MPI_Issend((void *) (write_send_buffer_size + proc_num), 1,
                      MPI_INT, proc_num, 0, world, write_size_speakers + proc_num);
            MPI_Issend((void *) (write_send_buffer + proc_num * BUFFER_SIZE), BUFFER_SIZE * sizeof(struct write_request),
                      MPI_BYTE, proc_num, 1, world, write_data_speakers + proc_num);
        }
        //omp_unset_lock(&write_send_locks[proc_num]);
    }
    return 1;
}

int empty_write_send_buffers() {
    int completed = 1;
    for (int proc_num = 0; proc_num < world_size; proc_num++) {
        if (proc_num == world_rank)
            continue;

        //omp_set_lock(&write_send_locks[proc_num]);
        if (write_send_buffer_size[proc_num] == BUFFER_SIZE) {
            int size_sent, data_sent;
            MPI_Test(write_size_speakers + proc_num, &size_sent, MPI_STATUS_IGNORE);
            MPI_Test(write_data_speakers + proc_num, &data_sent, MPI_STATUS_IGNORE);
            if (!size_sent || !data_sent) {
                completed = 0;
                //omp_unset_lock(&write_send_locks[proc_num]);
                continue;
            }
            write_send_buffer_size[proc_num] = 0;
        }

        if (write_send_buffer_size[proc_num] > 0) {
            MPI_Issend((void *) (write_send_buffer_size + proc_num),
                    1,
                    MPI_INT,
                    proc_num,
                    0,
                    world,
                    write_size_speakers + proc_num);
            MPI_Issend((void *) (write_send_buffer + proc_num * BUFFER_SIZE),
                      write_send_buffer_size[proc_num] * sizeof(struct write_request),
                              MPI_BYTE,
                              proc_num,
                              1,
                              world,
                              write_data_speakers + proc_num);
        }
        //omp_unset_lock(&write_send_locks[proc_num]);
    }
    if (completed) {
        close_write();
    }
    return completed;
}

void check_write_recv() {
    //omp_set_lock(&write_check_lock);
    if (world_rank > 0) {
        int count;
        int *indices = malloc(sizeof(int) * world_rank);
        MPI_Testsome(world_rank, write_listeners, &count, indices, MPI_STATUSES_IGNORE);
        for (int i = 0; i < count; i++) {
            if (open_active_write[indices[i]] >= 0)
                empty_write_recv_buffer(indices[i]);
        }
        free(indices);
    }
    if (world_rank < world_size - 1) {
        int count;
        int *indices = malloc(sizeof(int) * (world_size - world_rank - 1));
        MPI_Testsome(world_size - world_rank - 1, write_listeners + world_rank + 1, &count, indices, MPI_STATUSES_IGNORE);
        for (int i = 0; i < count; i++) {
            if (open_active_write[indices[i] + world_rank + 1] >= 0)
                empty_write_recv_buffer(indices[i] + world_rank + 1);
        }
        free(indices);
    }
    //omp_unset_lock(&write_check_lock);
}

void empty_write_recv_buffer(int proc_num) {
    if (write_recv_buffer_size[proc_num] == -1) {
        open_active_write[proc_num] = -1;
        return;
    } else if (write_recv_buffer_size[proc_num] == 0) {
        open_active_write[proc_num] = 0;
    } else {
        MPI_Recv((void *) (write_recv_buffer + proc_num * BUFFER_SIZE),
                 write_recv_buffer_size[proc_num] * sizeof(struct write_request), MPI_BYTE, proc_num, 1, world,
                 MPI_STATUS_IGNORE);
        for (int i = 0; i < write_recv_buffer_size[proc_num]; i++) {
            struct write_request request = write_recv_buffer[proc_num * BUFFER_SIZE + i];
            u64 h_f_x = murmur64(request.f_x);
            u64 indx = h_f_x % global_dict_size;

            u32 p_f_x = request.f_x % PRIME;

            u64 local_index = indx - min_local_h;
            //omp_set_lock(&write_recv_lock);
            while (local_dictionary[local_index].k != EMPTY) {
                local_index = (local_index + 1) % local_dict_size;
                if (local_index == indx - min_local_h) {
                    fprintf(stderr, "Error: no space in local dictionary");
                    exit(1);
                }
            }
            local_dictionary[local_index].k = p_f_x;
            local_dictionary[local_index].v = request.x;
            //omp_unset_lock(&write_recv_lock);
        }
        write_recv_buffer_size[proc_num] = 0;
    }
    MPI_Irecv(write_recv_buffer_size + proc_num, 1, MPI_INT, proc_num, 0, world, write_listeners + proc_num);
}

void empty_write_recv_buffers() {
    for (int i = 0; i < world_size; i++) {
        if (i == world_rank)
            continue;

        if (write_recv_buffer_size[i] > 0) {
            empty_write_recv_buffer(i);
        }
    }
}

int read_from_dict(u64 g_y, u64 y) {
    u64 h_g_y = murmur64(g_y);
    u64 indx = h_g_y % global_dict_size;
    int proc_num = indx / local_dict_size;

    u32 p_g_y = g_y % PRIME;

    if (proc_num == world_rank) {
        if (num_solutions == maxres)
            return 1;
        u64 local_index = indx - min_local_h;

        while (local_dictionary[local_index].k != EMPTY) {
            u64 p_f_x = local_dictionary[local_index].k;
            if (p_f_x == p_g_y) { // potential match
                u64 x = local_dictionary[local_index].v;
                local_candidates++;
                u64 f_x = f(x);
                if (f_x == g_y && is_good_pair(x,y)) { // solution
                    solutions[num_solutions].x = x;
                    solutions[num_solutions].y = y;
                    num_solutions++;
                    if (num_solutions == maxres)
                        break;
                }
            }
            local_index = (local_index + 1) % local_dict_size;
            if (local_index == indx - min_local_h) {
                break;
            }
        }
    } else {
        if (read_send_buffer_size[proc_num] == BUFFER_SIZE) {
            int size_sent, data_sent;
            MPI_Test(read_size_speakers + proc_num, &size_sent, MPI_STATUS_IGNORE);
            MPI_Test(read_data_speakers + proc_num, &data_sent, MPI_STATUS_IGNORE);
            if (!size_sent || !data_sent) {
                return 0;
            }
            read_send_buffer_size[proc_num] = 0;
        }
        struct read_request *addr = read_send_buffer + (proc_num * BUFFER_SIZE + read_send_buffer_size[proc_num]);
        addr->g_y = g_y;
        addr->y = y;
        read_send_buffer_size[proc_num]++;

        if (read_send_buffer_size[proc_num] == BUFFER_SIZE) {
            MPI_Issend((void *) (read_send_buffer_size + proc_num), 1,
                      MPI_INT, proc_num, 2, world, read_size_speakers + proc_num);
            MPI_Issend((void *) (read_send_buffer + proc_num * BUFFER_SIZE), BUFFER_SIZE * sizeof(struct read_request),
                      MPI_BYTE, proc_num, 3, world, read_data_speakers + proc_num);
        }
    }
    return 1;
}

int empty_read_send_buffers() {
    int completed = 1;
    for (int proc_num = 0; proc_num < world_size; proc_num++) {
        if (proc_num == world_rank)
            continue;

        if (read_send_buffer_size[proc_num] == BUFFER_SIZE) {
            int size_sent, data_sent;
            MPI_Test(read_size_speakers + proc_num, &size_sent, MPI_STATUS_IGNORE);
            MPI_Test(read_data_speakers + proc_num, &data_sent, MPI_STATUS_IGNORE);
            if (!size_sent || !data_sent) {
                completed = 0;
                continue;
            }
            read_send_buffer_size[proc_num] = 0;
        }

        if (read_send_buffer_size[proc_num] > 0) {
            MPI_Issend((void *) (read_send_buffer_size + proc_num),
                      1,
                      MPI_INT,
                      proc_num,
                      2,
                      world,
                      read_size_speakers + proc_num);
            MPI_Issend((void *) (read_send_buffer + proc_num * BUFFER_SIZE),
                      read_send_buffer_size[proc_num] * sizeof(struct read_request),
                      MPI_BYTE,
                      proc_num,
                      3,
                      world,
                      read_data_speakers + proc_num);
        }
    }
    if (completed) {
        close_read();
    }
    return completed;
}

void check_read_recv() {
    if (world_rank > 0) {
        int count;
        int *indices = malloc(sizeof(int) * world_rank);
        MPI_Testsome(world_rank, read_listeners, &count, indices, MPI_STATUSES_IGNORE);
        for (int i = 0; i < count; i++) {
            if (open_active_read[indices[i]] >= 0)
                empty_read_recv_buffer(indices[i]);
        }
        free(indices);
    }
    if (world_rank < world_size - 1) {
        int count;
        int *indices = malloc(sizeof(int) * (world_size - world_rank - 1));
        MPI_Testsome(world_size - world_rank - 1, read_listeners + world_rank + 1, &count, indices, MPI_STATUSES_IGNORE);
        for (int i = 0; i < count; i++) {
            if (open_active_read[indices[i] + world_rank + 1] >= 0)
                empty_read_recv_buffer(indices[i] + world_rank + 1);
        }
        free(indices);
    }
}

void empty_read_recv_buffer(int proc_num) {
    if (read_recv_buffer_size[proc_num] == -1) {
        open_active_read[proc_num] = -1;
        return;
    }
    if (read_recv_buffer_size[proc_num] == 0) {
        open_active_read[proc_num] = 0;
    } else {
        MPI_Recv((void *) (read_recv_buffer + proc_num * BUFFER_SIZE),
                 read_recv_buffer_size[proc_num] * sizeof(struct read_request), MPI_BYTE, proc_num, 3, world,
                 MPI_STATUS_IGNORE);
        //#pragma omp parallel for
        for (int i = 0; i < read_recv_buffer_size[proc_num]; i++) {
            if (num_solutions >= maxres)
                continue;
            struct read_request request = read_recv_buffer[proc_num * BUFFER_SIZE + i];
            u64 h_g_y = murmur64(request.g_y);
            u64 indx = h_g_y % global_dict_size;

            u32 p_g_y = request.g_y % PRIME;
            u64 local_index = indx - min_local_h;

            while (local_dictionary[local_index].k != EMPTY) {
                u64 p_f_x = local_dictionary[local_index].k;
                if (p_f_x == p_g_y) { // potential match
                    u64 x = local_dictionary[local_index].v;
                    //#pragma omp atomic
                    local_candidates++;
                    if (f(x) == request.g_y && is_good_pair(x, request.y) && num_solutions < maxres) { // solution
                        //#pragma omp critical
                        {
                            solutions[num_solutions].x = x;
                            solutions[num_solutions].y = request.y;
                            num_solutions++;
                        }
                    }
                }
                local_index = (local_index + 1) % local_dict_size;
                if (local_index == indx - min_local_h) {
                    break;
                }
            }
        }
        read_recv_buffer_size[proc_num] = 0;
    }
    MPI_Irecv(read_recv_buffer_size + proc_num, 1, MPI_INT, proc_num, 2, world, read_listeners + proc_num);
}

void empty_read_recv_buffers() {
    for (int i = 0; i < world_size; i++) {
        if (i == world_rank)
            continue;

        if (read_recv_buffer_size[i] > 0) {
            empty_read_recv_buffer(i);
        }
    }
}

void gather_solutions(int root) {
    int *nums;
    u64 global_candidates;
    struct result *sols;
    if (world_rank == root) {
        nums = malloc(sizeof(int) * world_size);
        sols = malloc(sizeof(struct result) * maxres * world_size);
    }
    MPI_Reduce((void *) &local_candidates, &global_candidates, 1, MPI_INT64_T, MPI_SUM, 0, world);
    MPI_Gather((void *) &num_solutions, 1, MPI_INT, nums, 1, MPI_INT, root, world);
    MPI_Gather((void *) solutions, 16 * sizeof(struct result), MPI_BYTE, sols, 16 * sizeof(struct result), MPI_BYTE, root, world);
    if (world_rank != root)
        return;
    int sum = 0;
    for (int i = 0; i < world_size; i++) {
        sum += nums[i];
    }
    printf(" %" PRId64 " candidate pairs tested\n", global_candidates);
    printf("Solutions found: %d\n", sum);
    for (int i = 0; i < world_size; i++) {
        for (int j = 0; j < nums[i]; j++) {
            printf("Solution found: (%" PRIx64 ", %" PRIx64 ") [checked OK] (f(%" PRIx64 ") = %" PRIx64 "; g(%" PRIx64 ") = %" PRIx64 ")\n",sols[i*maxres + j].x,sols[i*maxres + j].y,sols[i*maxres + j].x,f(sols[i*maxres + j].x),sols[i*maxres + j].y,g(sols[i*maxres + j].y));
        }
    }
    free(nums);
    free(sols);
}

void deactivate_write() {
    for (int i = 0; i < world_size; i++) {
        if (i == world_rank)
            continue;
        int num = 0;
        MPI_Send((void *) &num, 1, MPI_INT, i, 0, world);
    }
}

void deactivate_read() {
    for (int i = 0; i < world_size; i++) {
        if (i == world_rank)
            continue;
        int num = 0;
        MPI_Send((void *) &num, 1, MPI_INT, i, 2, world);
    }
}

void close_write() {
    for (int i = 0; i < world_size; i++) {
        if (i == world_rank)
            continue;
        int num = -1;
        MPI_Send((void *) (&num), 1, MPI_INT, i, 0, world);
    }
}

void close_read() {
    for (int i = 0; i < world_size; i++) {
        if (i == world_rank)
            continue;
        int num = -1;
        MPI_Send((void *) (&num), 1, MPI_INT, i, 2, world);
    }
}

int is_write_active() {
    for (int i = 0; i < world_size; i++) {
        if (open_active_write[i] == 1) {
            return 1;
        }
    }
    return 0;
}

int is_read_active() {
    for (int i = 0; i < world_size; i++) {
        if (open_active_read[i] == 1) {
            return 1;
        }
    }
    return 0;
}

int is_write_open() {
    for (int i = 0; i < world_size; i++) {
        if (open_active_write[i] >= 0) {
            return 1;
        }
    }
    return 0;
}

int is_read_open() {
    for (int i = 0; i < world_size; i++) {
        if (open_active_read[i] >= 0) {
            return 1;
        }
    }
    return 0;
}

void memory_cleanup() {
    free(local_dictionary);
    free(write_listeners);
    free(write_size_speakers);
    free(write_data_speakers);
    free(write_send_buffer_size);
    free(write_recv_buffer_size);
    free(write_send_buffer);
    free(write_recv_buffer);
    free(open_active_write);
    free(read_listeners);
    free(read_size_speakers);
    free(read_data_speakers);
    free(read_send_buffer_size);
    free(read_recv_buffer_size);
    free(read_send_buffer);
    free(read_recv_buffer);
    free(open_active_read);
    /*for (int i = 0; i < world_size; i++) {
        omp_destroy_lock(&write_send_locks[i]);
        omp_destroy_lock(&read_send_locks[i]);
    }
    omp_destroy_lock(&write_recv_lock);
    omp_destroy_lock(&write_check_lock);
    omp_destroy_lock(&read_recv_lock);
    omp_destroy_lock(&read_check_lock);
    free(write_send_locks);
    free(read_send_locks);*/
}

/* for debugging, please change buffer size for larger testing */
void print_state() {
    char arr[20000], *o = arr;
    sprintf(o,"%d>world rank: %d\n%d>global dict: %lu\n%d>min: %lu\n%d>max: %lu\n%d>local: %lu\n%d>dictionary: [\n", world_rank, world_rank, world_rank, global_dict_size, world_rank, min_local_h, world_rank, max_local_h, world_rank, local_dict_size, world_rank);
    o += strlen(o);
    for (u64 i = 0; i < local_dict_size; i++) {
        if (local_dictionary[i].k == EMPTY)
            sprintf(o, "%d>\t%lu (%lu): EMPTY\n", world_rank, i, i + min_local_h);
        else
            sprintf(o, "%d>\t%lu (%lu): {p_f_x: %d, x: %lu}\n", world_rank, i, i + min_local_h, local_dictionary[i].k, local_dictionary[i].v);
        o += strlen(o);
    }
    sprintf(o, "%d>]\n%d>write send buffer: [\n",world_rank,world_rank);
    o += strlen(o);
    for (int i = 0; i < world_size; i++) {
        sprintf(o, "%d>\t%d: [\n",world_rank,i);
        o += strlen(o);
        for (int j = 0; j < write_send_buffer_size[i]; j++) {
            sprintf(o, "%d>\t\t{f_x: %lu, x: %lu} (h(f_x) = %lu, h(f_x) %% global_dict_size = %lu)\n", world_rank, write_send_buffer[i*BUFFER_SIZE + j].f_x, write_send_buffer[i*BUFFER_SIZE + j].x,
                    murmur64(write_send_buffer[i*BUFFER_SIZE + j].f_x), murmur64(write_send_buffer[i*BUFFER_SIZE + j].f_x) % global_dict_size);
            o += strlen(o);
        }
        sprintf(o, "%d>\t]\n", world_rank);
        o += strlen(o);
    }
    sprintf(o, "%d>]\n%d>write recv buffer: [\n", world_rank, world_rank);
    o += strlen(o);
    for (int i = 0; i < world_size; i++) {
        sprintf(o, "%d>\t%d: [\n",world_rank,i);
        o += strlen(o);
        for (int j = 0; j < write_recv_buffer_size[i]; j++) {
            sprintf(o, "%d>\t\t{f_x: %lu, x: %lu} (h(f_x) = %lu, h(f_x) %% global_dict_size = %lu)\n", world_rank, write_recv_buffer[i*BUFFER_SIZE + j].f_x, write_recv_buffer[i*BUFFER_SIZE + j].x,
                    murmur64(write_recv_buffer[i*BUFFER_SIZE + j].f_x), murmur64(write_recv_buffer[i*BUFFER_SIZE + j].f_x) % global_dict_size);
            o += strlen(o);
        }
        sprintf(o, "%d>\t]\n", world_rank);
        o += strlen(o);
    }
    sprintf(o,"%d>]\n%d>read send buffer: [\n", world_rank, world_rank);
    o += strlen(o);
    for (int i = 0; i < world_size; i++) {
        sprintf(o, "%d>\t%d: [\n",world_rank,i);
        o += strlen(o);
        for (int j = 0; j < read_send_buffer_size[i]; j++) {
            sprintf(o, "%d>\t\t{g_y: %lu, y: %lu} (h(g_y) = %lu, h(g_y) %% global_dict_size = %lu)\n", world_rank, read_send_buffer[i*BUFFER_SIZE + j].g_y, read_send_buffer[i*BUFFER_SIZE + j].y,
                    murmur64(read_send_buffer[i*BUFFER_SIZE + j].g_y), murmur64(read_send_buffer[i*BUFFER_SIZE + j].g_y) % global_dict_size);
            o += strlen(o);
        }
        sprintf(o, "%d>\t]\n", world_rank);
        o += strlen(o);
    }
    sprintf(o,"%d>]\n%d>read recv buffer: [\n", world_rank, world_rank);
    o += strlen(o);
    for (int i = 0; i < world_size; i++) {
        sprintf(o, "%d>\t%d: [\n",world_rank,i);
        o += strlen(o);
        for (int j = 0; j < read_recv_buffer_size[i]; j++) {
            sprintf(o, "%d>\t\t{g_y: %lu, y: %lu} (h(g_y) = %lu, h(g_y) %% global_dict_size = %lu)\n", world_rank, read_recv_buffer[i*BUFFER_SIZE + j].g_y, read_recv_buffer[i*BUFFER_SIZE + j].y,
                    murmur64(read_recv_buffer[i*BUFFER_SIZE + j].g_y), murmur64(read_recv_buffer[i*BUFFER_SIZE + j].g_y) % global_dict_size);
            o += strlen(o);
        }
        sprintf(o, "%d>\t]\n", world_rank);
        o += strlen(o);
    }
    sprintf(o,"%d>]\n%d>solutions: [\n", world_rank, world_rank);
    o += strlen(o);
    for (int i = 0; i < num_solutions; i++) {
        sprintf(o, "%d>\t{x: %lu, y: %lu} (f(x) = %lu, g(x) = %lu)\n", world_rank, solutions[i].x, solutions[i].y, f(solutions[i].x), g(solutions[i].y));
        o += strlen(o);
    }
    sprintf(o,"%d>]\n", world_rank);
    o += strlen(o);

    printf("%s",arr);
}