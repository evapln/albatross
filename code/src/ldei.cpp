#include <cstdbool>
#include <NTL/ZZ_pX.h>
#include "ldei.hpp"

// in ZZq
bool localldei(const Vec<ZZ_p>& a, const long k, const Vec<ZZ_p>& x,
  const long m, const ZZ& p) {

  // computation of vector u
  Vec<ZZ_p> u;
  u.SetLength(m);
  ZZ_p prod, tmp;
  for (int i = 0; i < m; i++) {
    prod = 1;
    for (int l = 0; l < m; l++)
      if (l != i) {
        sub(tmp, a[i], a[l]);
        prod *= tmp;
      }
    inv(u[i],prod);
  }

  // random polynomial
  ZZ_pX P;
  do random(P,m-k-1); while(IsZero(P));

  // computation of v
  Vec<ZZ_p> v;
  v.SetLength(m);
  for (int i = 0; i < m; i++) {
    eval(tmp,P,a[i]);
    mul(v[i],u[i],tmp);
  }



  // verification
  ZZ_pPush push(p);
  prod = 1;
  for (int i = 0; i < m; i++) {
    power(tmp, x[i], rep(v[i]));
    prod *= tmp;
  }

  return (IsOne(prod));
}
