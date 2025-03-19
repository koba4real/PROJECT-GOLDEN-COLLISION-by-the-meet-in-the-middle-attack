#ifndef PROJET_UTIL_H
#define PROJET_UTIL_H

#include <stdbool.h>
#include <stdint.h>

typedef uint64_t u64;
typedef uint32_t u32;

extern u64 n;         /* block size (in bits) */
extern u64 mask;          /* this is 2**n - 1 */
/* (P, C) : two plaintext-ciphertext pairs */
extern u32 P[2][2];
extern u32 C[2][2];

void init_util();

void Speck64128KeySchedule(const u32 K[],u32 rk[]);

void Speck64128Encrypt(const u32 Pt[], u32 Ct[], const u32 rk[]);

void Speck64128Decrypt(u32 Pt[], const u32 Ct[], u32 const rk[]);

u64 murmur64(u64 x);

u64 f(u64 x);

u64 g(u64 y);

bool is_good_pair(u64 k1, u64 k2);

double wtime();

void process_command_line_options(int argc, char ** argv);

void usage(char **argv);

#endif //PROJET_UTIL_H
