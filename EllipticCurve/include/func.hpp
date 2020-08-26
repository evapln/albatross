#include <gmpxx.h>
extern "C" {
  #include "relic.h"
}

/* initiation of the randomness */
void prng_init(const unsigned int seed);

/* set y = P(x) with P polynomial of degree deg in Zq*/
void apply_poly(bn_t y, const bn_t *P, const int deg, const int x, bn_t q);
