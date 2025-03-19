#include <stdio.h>
#include <mpi.h>
#include <unistd.h>
#include "util.h"
#include "dictionary.h"

int main(int argc, char **argv) {
    MPI_Init(NULL,NULL);
    init_util();
    process_command_line_options(argc, argv);
    setup_dict(2 * (1ull << n), MPI_COMM_WORLD, 16);

    if (world_rank == 0) {
        printf("Running with n=%d, C0=(%08x, %08x) and C1=(%08x, %08x)\n",
               (int) n, C[0][0], C[0][1], C[1][0], C[1][1]);
        /*for (int i = 0; i < 16; i++) {
            printf("f(%d) = %lu\n", i, f(i));
        }
        for (int i = 0; i < 16; i++) {
            printf("g(%d) = %lu\n", i, g(i));
        }
        for (int i = 0; i < 16; i++) {
            for (int j = 0; j < 16; j++) {
                if (is_good_pair(i,j))
                    printf("%d, %d is good\n", i, j);
            }
        }*/
    } else sleep(1);

    //for (int i = world_rank * 4; i < 4 * (world_rank + 1); i++) {
    for (int i = world_rank; i < 16; i+=4) {
        write_to_dict(f(i),i);
    }
    empty_write_send_buffers();
    MPI_Barrier(MPI_COMM_WORLD);
    check_write_recv();
    MPI_Barrier(MPI_COMM_WORLD);

    //for (int i = world_rank * 4; i < 4 * (world_rank + 1); i++) {
    for (int i = world_rank; i < 16; i+=4) {
        read_from_dict(g(i),i);
    }
    empty_read_send_buffers();
    MPI_Barrier(MPI_COMM_WORLD);
    check_read_recv();
    //print_state();
    MPI_Barrier(MPI_COMM_WORLD);
    close_dict();
    gather_solutions(0);
    /*// manual 40 * 1.25
    setup_dict(50, MPI_COMM_WORLD, 16);
    for (int i = world_rank * 2; i < 2 * (world_rank + 1); i++) {
        write_to_dict(i+10,i);
    }
    sleep(world_rank + 1);
    print_state();
    sleep(5);

    empty_write_send_buffers();

    MPI_Barrier(world);

    sleep(world_rank + 1);

    check_write_recv();

    print_state();
*/
    MPI_Finalize();
}