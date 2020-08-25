#include <NTL/ZZ.h>

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
