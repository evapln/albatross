#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/matrix.h>

using namespace NTL;

/* initiation of the randomness */
void prng_init(const unsigned int seed);

/* set p and q such that:
  - q is prime
  - q > 2^(l+k)
  - 2^k divides q-1
  - p = 2p+1 is prime */
void findprime(ZZ& q, ZZ& p, const long k, const long l);

/* set w such that:
  - w is in ZZ_q
  - w is a n-th root of the unity */
void rootunity(ZZ& w, const long n, const ZZ& q);

/* compute a generator of the multiplicative group Zp* with p = 2q+1 */
void generator(ZZ_p& g, const ZZ& q);


////////////////////////////////////////////////////////////////////////////////
////////////////////// gestion of vectors and matrices /////////////////////////
////////////////////////////////////////////////////////////////////////////////

/* set the row i of M to vec */
int mat_set_row(Mat<ZZ_p>& M, const int i, const Vec<ZZ_p>& vec);

/* return true if vec is the null vector */
bool vec_is_zero(const Vec<ZZ_p>& vec);

/* exchange row1 and row2 in M */
void mat_exchange_row(Mat<ZZ_p>& M, const int row1, const int row2);

/* set row1 of M ot row1 + coef * row2 */
void mat_add_row(Mat<ZZ_p>& M, const int row1, const int row2, const ZZ_p& coef);

/* systematise M  on its 'row' first rows */
void mat_systematisation(Mat<ZZ_p>& M, const int row);

/* return true if the 'row' first rows of M are independants, false otherwise */
bool is_ind(const Mat<ZZ_p>& M, const int row);
