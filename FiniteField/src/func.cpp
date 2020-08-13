#include <cstdlib>
#include <stdbool.h>

#include "func.hpp"


void prng_init(const unsigned int seed) {
  static bool seed_init = false;
  if (! seed_init) {
    srand(seed);
    seed_init = true;
  }
}

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
