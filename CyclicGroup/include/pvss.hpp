#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/vector.h>
#include <NTL/matrix.h>

using namespace std;
using namespace NTL;

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

/* compute the list of Lagrange coefficients */
void lambda(Mat<ZZ_p>& lamb, const int t, pl_t *pl);

/* create the public ledger, create the public and private keys for n participants
return the time used to compute the exponentiations */
clock_t setup(pl_t* pl, Vec<ZZ_p>& sk, const int n, const ZZ& q, const ZZ& p, const ZZ_p& h);

/* add the encrypted shares in he public ledger with t threshold, l number of secrets
return the time used to compute the exponentiations */
clock_t distribution(const int l, const int t, const Vec<ZZ_p>& alpha, pl_t *pl);

/* reconstruct the secret vector
return the time used to compute the exponentiations */
clock_t reconstruction(const int r, pl_t *pl);


////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// test /////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void pvss_test(const int n, const int size);
