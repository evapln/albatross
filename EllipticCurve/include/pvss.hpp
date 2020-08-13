#include <gmpxx.h>
extern "C" {
  #include "relic.h"
}

using namespace std;

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// public ledger /////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

typedef struct pl_t pl_t;

/* allocate the memory for the publice ledger */
pl_t *pl_alloc(const int n);

/* liberate the memory of the public ledger */
void pl_free(pl_t *pl);

/* print the public ledger */
void pl_print(pl_t *pl);


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// scheme //////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/* setup return the public ledger
and create the public keys and secret keys for n participants */
pl_t *setup(bn_t* sk, const int n, bn_t q);

/* set y = P(x) with P polynomila of degree deg */
void apply_poly(bn_t y, const bn_t *P, const int deg, const int x, bn_t q);

/* t threshold, l number of secrets
add the encrypted shares in he public ledger */
void distribution(const int l, const int t, pl_t *pl, bn_t q);

/* compute the list of lambda to compute the secrets */
void lambda(bn_t** lambs, const int t, pl_t *pl, bn_t q);

/* reconstruct the secret vector */
void reconstruction(const int r, pl_t *pl, bn_t q);

void pvss_test(void);
