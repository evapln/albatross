#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/vector.h>

using namespace std;
using namespace NTL;

typedef struct pl_t pl_t;

void prng_init(const unsigned int seed);

pl_t *pl_alloc(const int n);
void pl_free(pl_t *pl);
void pl_print(pl_t *pl);

void generator(ZZ_p& g, const ZZ& q);

// setup return the secretkey
// float *setup(const int q, const double h, pl_t *pl);
pl_t *setup(Vec<ZZ_p>& sk, const int n, const ZZ& q, const ZZ& p, const ZZ_p& h);

// n participants, t threshold, l = n-2t-1
// return the list of the secrets ( one secret for one participant)
void distribution(const int l, const int t, pl_t *pl);

// int lambda(const int i, const int j, const int *tab, const int len);
void lambda(Mat<ZZ_p>& lamb, int t, pl_t *pl);

// reconstruct the secret vector (well... now no... but soon...)
void reconstruction(Vec<ZZ>& S, const int r, pl_t *pl);
