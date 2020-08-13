/* Implementation of albatross over tweedledum elliptic curve defined as:
  Ep : y^2 = x^3 + 5 over GF(p) of order q with 
    p = 2^254 + 4707489545178046908921067385359695873
    q = 2^254 + 4707489544292117082687961190295928833
*/

#include "pvss.hpp"

#include <iostream>
#include <cstdlib>
#include <ctime>

using namespace std;

int main(void) {
  cout << "pvss" << endl;
  pvss_test();
  return EXIT_SUCCESS;
}
