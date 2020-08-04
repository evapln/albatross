#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/vector.h>
#include <NTL/matrix.h>


using namespace NTL;

// (q, p = 2*q+1) such that q is prime, q>2^(l+k), 2*q+1 is prime and 2^k divides q-1
void findprime(const long k, const long l, ZZ& q, ZZ& p);
void rootunity(ZZ& w, const long n, const ZZ& q);
void FFTE(Vec<ZZ_p>& f, const long n, const Vec<ZZ_p>& h, const ZZ& w, const ZZ& q);
void evalu(ZZ& sum, const Vec<ZZ_p>& L, const ZZ& x, const ZZ& q);
bool test(const Vec<ZZ_p>& W, const Vec<ZZ_p>& coef, const ZZ_p& h, const ZZ& w,
          const long n, const ZZ& q);


int mat_set_row(Mat<bool>& M, const int i, const Vec<bool>& vec);
bool vec_is_zero(const Vec<ZZ_p>& vec);
// bool vec_is_in_mat(const Mat<bool>& M, const Vec<bool>& v);

// int space_add_row(Mat<bool>& S, const int i, const Vec<bool>& r);

void mat_exchange_row(Mat<ZZ_p>& M, const int row1, const int row2);
void mat_sub_row(Mat<ZZ_p>& M, const int row1, const int row2);
void mat_systematisation(Mat<ZZ_p>& M, const int row);
bool is_ind(const Mat<ZZ_p>& M, const int row);
// return a non zero binary vector of length n
Vec<ZZ_p> bin_rand_vec(const int n);
void bin_rand_vec(Vec<ZZ_p>& vec, const int n);

// int mat_gen_space(Mat<bool>& M, const int n, const int k);
void mat_gen(Mat<ZZ_p>& M, const int n, const int k);

void mul_mat(Vec<ZZ_p>& hhat, Mat<ZZ_p>& M, const int n, const int k, const Vec<ZZ_p>& h);
