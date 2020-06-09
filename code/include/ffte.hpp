#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/vector.h>


using namespace NTL;

// (q, p = 2*q+1) such that q is prime, q>2^(l+k), 2*q+1 is prime and 2^k divides q-1
void findprime(const long k, const long l, ZZ& q, ZZ& p);
void rootunity(ZZ& w, const long n, const ZZ& q);
void FFTE(Vec<ZZ_p>& f, const long n, const Vec<ZZ_p>& h, const ZZ& w, const ZZ& q);
void evalu(ZZ& sum, const Vec<ZZ_p>& L, const ZZ& x, const ZZ& q);
bool test(const Vec<ZZ_p>& W, const Vec<ZZ_p>& coef, const ZZ_p& h, const ZZ& w,
          const long n, const ZZ& q);

/* time for one call to FFTE with q of size 1024 and n=1024 -> about 2s */
