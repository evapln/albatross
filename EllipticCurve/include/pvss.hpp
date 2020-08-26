#include <gmpxx.h>
extern "C" {
  #include "relic.h"
}

using namespace std;

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// public ledger /////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

typedef struct pl_t pl_t;

/* allocate the memory for the public ledger */
pl_t *pl_alloc(const int n);

/* liberate the memory of the public ledger */
void pl_free(pl_t *pl);

/* print the public ledger */
void pl_print(pl_t *pl);


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// scheme //////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/* compute the list of Lagrange coefficients */
void lambda(bn_t** lambs, const int t, pl_t *pl, bn_t q);

/* create the public ledger, create the public and private keys for n participants
return the time used to compute the exponentiations */
clock_t setup(pl_t* pl, bn_t* sk, const int n, bn_t q);

/* add the encrypted shares in he public ledger with t threshold, l number of secrets
return the time used to compute the exponentiations */
clock_t distribution(const int l, const int t, pl_t *pl, bn_t q);

/* reconstruct the secret vector
return the time used to compute the exponentiations */
clock_t reconstruction(const int r, pl_t *pl, bn_t q);

////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// ppvss test ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void pvss_test(const int n);
