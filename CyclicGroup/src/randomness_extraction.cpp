#include <ctime>
#include <iostream>
#include <NTL/ZZ.h>
#include <math.h>
#include <stdbool.h>
#include <unistd.h>

#include "randomness_extraction.hpp"
#include "func.hpp"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
///////////////////// extraction of randomness: FFTE ///////////////////////////
////////////////////////////////////////////////////////////////////////////////

void FFTE(Vec<ZZ_p>& f, const long n, const Vec<ZZ_p>& h, const ZZ& w, const ZZ& q) {
  f.SetLength(n);
  if (n == 1) {
    f[0] = h[0];
    return;
  }
  long m = n/2;
  Vec<ZZ_p> u;
  u.SetLength(m);
  Vec<ZZ_p> v;
  v.SetLength(m);
  ZZ_p hj, hjm, invhjm, tmp;
  ZZ wj, w2;
  // computation of the vectors u and v
  for (int j = 0; j < m; j++) {
    hj = h[j];
    hjm = h[j+m];
    u[j] = hj*hjm;
    inv(invhjm,hjm);
    PowerMod(wj,w,j,q);
    mul (tmp, hj, invhjm);
    power(v[j], tmp, wj);
  }
  Vec<ZZ_p> uhat;
  uhat.SetLength(m);
  Vec<ZZ_p> vhat;
  vhat.SetLength(m);
  PowerMod(w2,w,2,q);
  // recursive calls FFTE(n/2,u,w^2) FFTE(n/2,v,w^2)
  FFTE(uhat,m,u,w2,q);
  FFTE(vhat,m,v,w2,q);
  // set of f
  for (int i = 0; i < m; i++) {
    f[i*2] = uhat[i];
    f[i*2+1] = vhat[i];
  }
}

void evalu(ZZ& sum, const Vec<ZZ_p>& L, const ZZ& x, const ZZ& q) {
  sum = 0;
  ZZ tmp, tmp1;
  for (int i = 0; i < L.length(); i++) {
    PowerMod(tmp, x, i, q);
    MulMod(tmp1, rep(L[i]), tmp, q);
    AddMod(sum, sum, tmp1, q);
  }
}

bool test(const Vec<ZZ_p>& W, const Vec<ZZ_p>& coef, const ZZ_p& h, const ZZ& w,
const long n, const ZZ& q) {
  ZZ eval;
  ZZ wi;
  Vec<ZZ_p> U;
  U.SetLength(n);
  for (int i = 0; i < n; i++) {
    PowerMod(wi,w,i,q);
    evalu(eval, coef, wi, q);
    power(U[i],h, eval);
  }
  return (U == W);
}


////////////////////////////////////////////////////////////////////////////////
/////////////////// extraction of randomness: second method ////////////////////
////////////////////////////////////////////////////////////////////////////////

// Vec<ZZ_p> bin_rand_vec(const int n) {
//   Vec<ZZ_p> vec;
//   vec.SetLength(n);
//   prng_init(time(NULL) + getpid());
//   do {
//     for (int i = 0; i < n; i++)
//       vec[i] = ZZ_p(rand() % 2);
//   } while (vec_is_zero(vec));
//   return vec;
// }

void bin_rand_vec(Vec<ZZ_p>& vec, const int n) {
  vec.SetLength(n);
  prng_init(time(NULL) + getpid());
  do {
    for (int i = 0; i < n; i++)
      vec[i] = ZZ_p(rand() % 2);
  } while (vec_is_zero(vec));
}

void mat_gen(Mat<ZZ_p>& M, const int n, const int k) {
  M.SetDims(k,n);
  Vec<ZZ_p> r;
  bool ind = false;
  while (!ind) {
    for (int i = 0; i < k; i++) {
      bin_rand_vec(r,n);
      mat_set_row(M,i,r);
    }
    ind = is_ind(M,k);
  }
}

void mul_mat(Vec<ZZ_p>& hhat, Mat<ZZ_p>& M, const int n, const int k, const Vec<ZZ_p>& h) {
  if (h.length() != n) {
    hhat.SetLength(0);
    return;
  }
  hhat.SetLength(k);
  for (int i = 0; i < k; i++) {
    hhat[i] = ZZ_p(1);
    for (int j = 0; j < n; j++) {
      if (M[i][j] == 1)
        hhat[i] *= h[j];
    }
  }
}
