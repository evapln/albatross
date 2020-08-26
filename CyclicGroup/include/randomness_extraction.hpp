#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/vector.h>
#include <NTL/matrix.h>

using namespace NTL;

////////////////////////////////////////////////////////////////////////////////
///////////////////// extraction of randomness: FFTE ///////////////////////////
////////////////////////////////////////////////////////////////////////////////

/* set the vector f as FFTE(n,h,w) with:
  - n = 2^k | q-1
  - h vector of Gq
  - w n-th root of the unity
  - q large prime order of Gq */
void FFTE(Vec<ZZ_p>& f, const long n, const Vec<ZZ_p>& h, const ZZ& w, const ZZ& q);

/* set sum at L(x), L is a polynomial modulo q */
void evalu(ZZ& sum, const Vec<ZZ_p>& L, const ZZ& x, const ZZ& q);

/* return true if W is equal to ffte(n,h,w), false otherwise */
bool test(const Vec<ZZ_p>& W, const Vec<ZZ_p>& coef, const ZZ_p& h, const ZZ& w, const long n, const ZZ& q);

////////////////////////////////////////////////////////////////////////////////
/////////////////// extraction of randomness: second method ////////////////////
////////////////////////////////////////////////////////////////////////////////

/* set vec to a non zero binary vector of Zp of length n */
void bin_rand_vec(Vec<ZZ_p>& vec, const int n);

/* set M to a random matrix of size nxk on Zp where all rows are independant */
void mat_gen(Mat<ZZ_p>& M, const int n, const int k);

/* set hhat at M \diamond h ie hhat_i = prod_{k=1}^n h_k^M_{ik}*/
void mul_mat(Vec<ZZ_p>& hhat, Mat<ZZ_p>& M, const int n, const int k, const Vec<ZZ_p>& h);
