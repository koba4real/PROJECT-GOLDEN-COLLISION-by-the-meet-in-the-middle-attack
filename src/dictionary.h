#ifndef PROJET_DICTIONARY_H
#define PROJET_DICTIONARY_H

#include <stdint.h>
#include <stdio.h>
#include <mpi.h>
#include <omp.h>
#include "util.h"

#define BUFFER_SIZE 1000
#define MAX_RESULTS 256
static const u32 EMPTY = 0xffffffff;
static const u64 PRIME = 0xfffffffb;

extern int maxres;

extern MPI_Comm world;
extern int world_size;
extern int world_rank;
extern u64 global_dict_size;
extern u64 min_local_h, max_local_h;
extern u64 local_dict_size;
struct __attribute__ ((packed)) entry { u32 k; u64 v; };  /* hash table entry */
extern struct entry *local_dictionary; // array

struct write_request { u64 f_x; u64 x; }; /* to send to other process */

extern MPI_Request *write_listeners;
extern MPI_Request *write_size_speakers;
extern MPI_Request *write_data_speakers;

extern struct write_request *write_send_buffer;
extern int *write_send_buffer_size;

extern struct write_request *write_recv_buffer;
extern int *write_recv_buffer_size;

extern MPI_Request *read_listeners;
extern MPI_Request *read_size_speakers;
extern MPI_Request *read_data_speakers;

struct read_request { u64 g_y; u64 y; };

extern struct read_request *read_send_buffer;
extern int *read_send_buffer_size;

extern struct read_request *read_recv_buffer;
extern int *read_recv_buffer_size;

struct result { u64 x; u64 y; }; // first entry contains query's x,
// second entry contains count of results (0 - MAX_RESULTS, or MAX_RESULTS + 1 for overflow),
// next MAX_RESULTS entries are results or padding

extern struct result solutions[MAX_RESULTS];
extern int num_solutions;
extern u64 local_candidates;

extern int *open_active_write; // -1 = closed, 0 = inactive, 1 = open
extern int *open_active_read;

/*extern omp_lock_t *write_send_locks;
extern omp_lock_t write_recv_lock;
extern omp_lock_t write_check_lock;
extern omp_lock_t *read_send_locks;
extern omp_lock_t read_recv_lock;
extern omp_lock_t read_check_lock;
extern omp_lock_t sol_lock;*/

void setup_dict(u64 global_dict_size, MPI_Comm comm_world, int max_results);

int write_to_dict(u64 f_x, u64 x);
int empty_write_send_buffers();
void check_write_recv();
void empty_write_recv_buffer(int proc_num);
void empty_write_recv_buffers();

int read_from_dict(u64 g_y, u64 y);
int empty_read_send_buffers();
void check_read_recv();
void empty_read_recv_buffer(int proc_num);
void empty_read_recv_buffers();

void gather_solutions(int rank);
void close_write();
void close_read();
int is_write_active();
int is_read_active();
void deactivate_write();
void deactivate_read();
int is_write_open();
int is_read_open();
void print_state();

void memory_cleanup();

#endif //PROJET_DICTIONARY_H
