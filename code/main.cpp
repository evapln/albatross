#include "ffte.h"
#include "pvss.h"

int main(void) {
  long n = 128;
  long k = 7;
  long l = 1024 - k;
  ZZ p,q;
  findprime(k,l,q,p);
  cout << "q = " << q << "  p = " << p << endl;
  ZZ w;
  rootunity(w,n,q);
  cout << "w = " << w << endl;
  ZZ_p::init(p);
  ZZ_p h {4};
  Vec<ZZ_p> L, coef;
  L.SetLength(n);
  coef.SetLength(n);
  for (int i = 0; i < n; i++) {
    random(coef[i]);
    power(L[i],h,rep(coef[i]));
  }
  Vec<ZZ_p> f;
  FFTE(f,n,L,w,q);
  cout << "f = " << f << endl;
  cout << "correct ? " << test(f,coef,h,w,n,q) << endl;
  return EXIT_SUCCESS;
}

/* TODO :
-> find the error in ffte
*/
