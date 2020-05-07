#ifndef PVSS_H
#define PVSS_H

typedef struct pl_t pl_t;

pl_t *pl_alloc(int n);
void pl_free(pl_t *pl);

int mod(int n, int q);

void prng_init(unsigned int seed);

int apply_poly(int *poly, int deg, int val, int mod);
int *rand_poly(int deg, int mod);

int *setup(int q, int h, pl_t *pl);
void distribution(int l, int t, pl_t *pl);
int lambda(int i, int j, int *tab, int len);
int *reconstruction(int *parts, int len_parts, pl_t *pl);

#endif
