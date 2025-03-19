#include <inttypes.h>
#include <stdio.h>
#include <mpi.h>
#include <unistd.h>
#include <omp.h>
#include "util.h"
#include "dictionary.h"

int main(int argc, char **argv) {
    double start = 0, mid = 0;

    MPI_Init(NULL,NULL);
    init_util();
    process_command_line_options(argc, argv);
    u64 N = 1ull << n;
    setup_dict(1.125 * N, MPI_COMM_WORLD, 16);

    if (world_rank == 0) {
        printf("Running with n=%d, C0=(%08x, %08x) and C1=(%08x, %08x)\n",
               (int) n, C[0][0], C[0][1], C[1][0], C[1][1]);
        start = wtime();
    }
    // compute
    for (u64 x = world_rank; x < N; x += world_size) {
        if ((x / world_size) % (1000) == 0) {
            check_write_recv();
        }
        while (!write_to_dict(f(x),x)) {
            check_write_recv();
        }
    }
    // signal computation finished
    deactivate_write();
    // empty buffer until all computation finished (otherwise can be blocking)
    while (is_write_active()) {
        check_write_recv();
    }
    // empty remains of buffer (and close)
    while(1) {
        if (empty_write_send_buffers())
            break;
    }
    // process buffers
    while (is_write_open()) {
        check_write_recv();
    }
    if (world_rank == 0) {
        mid = wtime();
        printf("Fill: %.1fs\n", mid - start);
    }
    // compute g(y) for checking solutions
    for (u64 y = world_rank; y < N; y += world_size) {
        if ((y / world_size) % (1000) == 0) {
            check_read_recv();
        }
        while (!read_from_dict(g(y),y)) {
            check_read_recv();
        }
    }
    // signal end of computation
    deactivate_read();
    // keep emptying buffers until computation finished (avoid deadlocks)
    while (is_read_active()) {
        check_read_recv();
    }
    // empty last of buffer
    while(1) {
        if (empty_read_send_buffers())
            break;
    }
    // receive last of buffer
    while (is_read_open()) {
        check_read_recv();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    // finished
    if (world_rank == 0) {
        printf("Probe: %.1fs.", wtime() - mid);
    }
    gather_solutions(0);
    memory_cleanup();
    MPI_Finalize();
}