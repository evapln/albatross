#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/vector.h>

using namespace std;
using namespace NTL;

typedef struct pl_t pl_t;

pl_t *pl_alloc(const int n);
void pl_free(pl_t *pl);
void pl_print(pl_t *pl);

int mod(const int n, const int q);

// initialisation of randomness -> call with seed = time(NULL) + getpid()
void prng_init(const unsigned int seed);

double apply_poly(const int *poly, const int deg, const int val, const int q);
// random polynomial in Z_q[X] of degree deg
int *rand_poly(const int deg, const int q);

// setup return the secretkey
float *setup(const int q, const double h, pl_t *pl);

// n participants, t threshold, l = n-2t-1
// return the list of the secrets ( one secret for one participant)
void distribution(const int l, const int t, pl_t *pl);
int lambda(const int i, const int j, const int *tab, const int len);

// reconstruct the secret vector (well... now no... but soon...)
void reconstruction(Vec<ZZ>& s, const int r, pl_t *pl);
