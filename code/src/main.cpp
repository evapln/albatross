#include "../include/ffte.hpp"
#include "../include/pvss.hpp"
#include "../include/ldei.hpp"


#include <iostream>
#include <cstdlib>
#include <ctime>
#include <NTL/ZZ_pX.h>

bool IsIn(const Vec<ZZ_p>& v, const int n, const ZZ_p val) {
    for (int i = 0; i < n; i++)
      if (v[i] == val)
        return true;
    return false;
}

Vec<ZZ_p> vect_rand_dist(const int n) {
  Vec<ZZ_p> v;
  v.SetLength(n);
  ZZ_p tmp;
  int len = 0;
  while (len < n) {
    random(tmp);
    if (!IsIn(v, len, tmp)) {
      v[len] = tmp;
      len++;
    }
  }
  return v;
}

int main(void) {
  // // test ffte
  // long n = 1024;
  // long k = 102;
  // long l = 1024 - k;
  // ZZ p,q;
  // findprime(k,l,q,p);
  // cout << "q = " << q << "  p = " << p << endl;
  // ZZ w;
  // rootunity(w,n,q);
  // cout << "w = " << w << endl;
  // ZZ_p::init(p);
  // ZZ_p h {4};
  // Vec<ZZ_p> L, coef;
  // L.SetLength(n);
  // coef.SetLength(n);
  // for (int i = 0; i < n; i++) {
  //   random(coef[i]);
  //   power(L[i],h,rep(coef[i]));
  // }
  // Vec<ZZ_p> f;
  // clock_t FFTE_time = clock();
  // FFTE(f,n,L,w,q);
  // FFTE_time = clock() - FFTE_time;
  // cout << "f = " << f << endl;
  // cout << "correct ? " << test(f,coef,h,w,n,q) << endl;
  // cout << "time: " << (float)FFTE_time/CLOCKS_PER_SEC << "s" << endl;

  // // test pvss
  // pvss();

  // test LDEI
  long m = 1024, k = 764;
  ZZ q = GenGermainPrime_ZZ(1024);
  ZZ p = 2 * q + 1;
  cout << "q = " << q << "  p = " << p << endl;
  ZZ_p::init(q);
  Vec<ZZ_p> a = vect_rand_dist(m);
  ZZ_pX P;
  random(P,k+1);
  Vec<ZZ_p> pi = eval(P,a);
  ZZ_p g;
  generator(g,q);

  ZZ_p::init(p);
  ZZ_p h = power(g,2);
  Vec<ZZ_p> x;
  x.SetLength(m);
  for (int i = 0; i < m; i++)
    power(x[i],h,rep(pi[i]));

  Vec<ZZ_p> xf;
  xf.SetLength(m);
  for (int i = 0; i < m; i++)
    random(xf[i]);

  ZZ_p::init(q);
  clock_t rec, time1, time2;
  bool verif1, verif2;

  rec = clock();
  verif1 = localldei(a, k, x, m, p);
  time1 = clock() - rec;

  rec = clock();
  verif2 = localldei(a, k, xf, m, p);
  time2 = clock() - rec;

  cout << "It's " << verif1 << " in " << (float)time1/CLOCKS_PER_SEC << " second(s).\n";
  cout << "It's " << verif2 << " in " << (float)time2/CLOCKS_PER_SEC << " second(s).\n";

  // return
  return EXIT_SUCCESS;
}
