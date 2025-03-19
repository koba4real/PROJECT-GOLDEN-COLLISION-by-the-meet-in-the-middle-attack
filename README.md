# PROJECT-GOLDEN-COLLISION-by-the-meet-in-the-middle-attack
project i develloped as part of my master degree in parallel programming.

The purpose of this project was to implement a parallelized version of a meet-in-the-middle attack program to find a "golden collision" for the problem defined as follows:

  Given two unary functions f, g mapping an n-bit integer to n-bit integers (assumed to be cryptographic hash functions, i.e. quick evaluation, internal analysis infeasible), and a binary relation π on pairs of n-bit integers, a golden collision is defined as a pair of n-bit integers (x, y) such that f(x) = g(y) and π(x, y) = 1. 
The sequential program provided (sequential.c), computes a brute force inverse of f in the form of a dictionary, with which g can be used, for a given y, to find values x such that f(x) = g(y), after which π(x,y) can be verifi ed. This serves as the baseline for parallelisation. This complexity of this approach is exponential, requiring at a minimum to compute 2^n evaluations of f, and then up to 2^n evaluations of g. The memory complexity is again exponential, requiring the storage of 2^n inputs to f in the dictionary.

The primary objective of parallelisation was to decrease the computation time and memory required per core, to increase the maximum input size, while also reducing total execution time.
