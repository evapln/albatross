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

void apply_poly(bn_t y, const bn_t *P, const int deg, const int x, bn_t q) {
  bn_t bnx, pow, tmp;
  bn_new(bnx);
  bn_new(y);
  bn_new(pow);
  bn_null(tmp);
  bn_copy(y, P[0]);
  if (x < 0) {
    bn_set_dig(bnx, -x);
    bn_neg(bnx,bnx);
  } else
    bn_set_dig(bnx, x);
  bn_copy(pow,bnx);
  for (int i = 1; i < deg; ++i) {
    bn_new(tmp);
    bn_mul(tmp, P[i], pow);
    bn_add(y,y,tmp);
    bn_mod(y,y,q);
    bn_mul(pow,pow,bnx);
    bn_mod(pow,pow,q);
  }
  bn_free(fpx);
  bn_free(pow);
  bn_free(tmp);
}
