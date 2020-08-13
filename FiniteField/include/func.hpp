#include <NTL/ZZ.h>

using namespace NTL;

void prng_init(const unsigned int seed);

// (q, p = 2*q+1) such that q is prime, q>2^(l+k), 2*q+1 is prime and 2^k divides q-1
void findprime(const long k, const long l, ZZ& q, ZZ& p);

void rootunity(ZZ& w, const long n, const ZZ& q);
