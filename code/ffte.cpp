#include <iostream>
#include <NTL/ZZ.h>
#include <math.h>
#include <stdbool.h>
#include "ffte.h"

using namespace std;

void findprime(const long k, const long l, ZZ& q, ZZ& p) {
  ZZ n, tmp;
  power2(n, k);
  long s = k % 2 - l % 2;
  power2(tmp, l); // tmp = 2**l
  q = (tmp + s)* n + 1; // q = n*(2**l+s) + 1
  p = 2 * q + 1;
  int r = 0;
  int bound = pow(10,8);
  while (r < bound) {
    if (ProbPrime(q) == 1 && ProbPrime(p) == 1)
      return;
    q += 3 * n;
    p += 6 * n;
    r++;
  }
}

// mod q
void rootunity(ZZ& w, const long n, const ZZ& q) {
  int i = 2;
  ZZ t = (q-1) / n;
  ZZ tmp;
  while (true) {
    PowerMod(w,ZZ(i),t,q);
    PowerMod(tmp, w, n/2,q);
    if (!IsOne(tmp))
      return;
    i++;
  }
}

// dans ZZp
void FFTE(Vec<ZZ_p>& f, const long n, const Vec<ZZ_p>& h, const ZZ& w, const ZZ& q) {
  f.SetLength(n);
  if (n == 1) {
    f[0] = h[0];
    return;
  }
  long m = n/2;
  Vec<ZZ_p> u; //vec de ZZp
  u.SetLength(m);
  Vec<ZZ_p> v; // vec de ZZp
  v.SetLength(m);
  ZZ_p hj, hjm, invhjm, tmp;
  ZZ wj, w2; // in ZZq
  for (int j = 0; j < m; j++) {
    hj = h[j];
    hjm = h[j+m];
    u[j] = hj*hjm;
    inv(invhjm,hjm);
    PowerMod(wj,w,j-1,q);
    mul (tmp, hj, invhjm);
    power(v[j], tmp, wj);
  }
  Vec<ZZ_p> uhat;
  uhat.SetLength(m);
  Vec<ZZ_p> vhat;
  vhat.SetLength(m);
  PowerMod(w2,w,2,q);
  FFTE(uhat,m,u,w2,q);
  FFTE(vhat,m,v,w2,q);
  for (int i = 0; i < m; i++) {
    f[i*2] = uhat[i];
    f[i*2+1] = vhat[i];
  }
}

// dans ZZq
void evalu(ZZ& sum, const Vec<ZZ_p>& L, const ZZ& x, const ZZ& q) {
  sum = 0;
  ZZ tmp, tmp1;
  for (int i = 0; i < L.length(); i++) {
    PowerMod(tmp, x, i, q);
    MulMod(tmp1, rep(L[i]), tmp, q);
    AddMod(sum, sum, tmp1, q);
  }
}

// dans ZZp
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
