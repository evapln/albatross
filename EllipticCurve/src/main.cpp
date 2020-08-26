/* Implementation of ppvss over tweedledum elliptic curve defined as:
  Ep : y^2 = x^3 + 5 over GF(p) of order q with
    p = 2^254 + 4707489545178046908921067385359695873
    q = 2^254 + 4707489544292117082687961190295928833
*/

#include "pvss.hpp"

#include <iostream>
#include <cstdlib>
#include <ctime>

using namespace std;

int main(int argc, char *argv[]) {
  int n = 200;
  if (argc == 2)
    n = atoi(argv[1]);
  cout << "ppvss test for " << n << " participants\n\n";
  pvss_test(n);
  return EXIT_SUCCESS;
}
