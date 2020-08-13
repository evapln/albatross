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
