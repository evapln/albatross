extern "C" {
  #include "relic/relic.h"
}
#include <fstream>
#include <iostream>

using namespace std;


#define TDUM_A		"0"
#define TDUM_B		"5"
#define TDUM_X		"40000000000000000000000000000000038AA1276C3F59B9A14064E200000000"
#define TDUM_Y		"2"
#define TDUM_R
#define TDUM_H
#define TDUM_MAPU	"5"

int power[44] = {254,121,120,119,115,113,111,109,104,101,98,97,96,94,93,91,90,85,84,83,82,81,80,78,76,75,72,71,69,68,67,64,63,61,56,54,46,45,42,39,38,37,33,0};

void bn_set_2_power(bn_t p, const int *f, const int len) {
  bn_t t;

	bn_null(t);
	bn_new(t);

  if (len < 1) {
    p = 0;
    return;
  }

	bn_set_2b(p, f[len - 1]);
	for (int i = len - 2; i > 0; i--) {
		bn_set_2b(t, f[i]);
		bn_add(p, p, t);
	}
}

int main(void) {
  // // opening of p.txt
  // ifstream pfile;
  // pfile.open("p.txt");
  // if (!pfile) {
  //   cout << "I'm not able to open the file\n";
  //   return EXIT_FAILURE;
  // }
  // // declaration of Fp
  // string pstr;
  // pfile >> pstr;
  // pfile.close();
  // int len = pstr.length();
  // const char *pchar = pstr.c_str();
  // cout << pstr << endl;
  // bn_t p;
  // bn_null(p);
  // bn_new(p);
  // bn_read_str(p,pchar,len,10);
  // pfile >> p;
  // p << pfile;
  // char pp[len];
  // bn_write_str(pp, len, p, 10);
  // cout << p << endl;
  // fp_prime_init();
  // fp_prime_set_dense(p);
  //
  // //declaration of the elliptic curve
  // fp_t a,b,u;
  // bn_t r,h;
  // ep_t g;
  // int ctmap = false;
  // ep_curve_set_plain(a, b, g, r, h, u, ctmap);
  //
  if (core_init() != RLC_OK) {
		core_clean();
		return 1;
	}


  // bn_t p;
  // bn_null(p);
  // bn_new(p);
  // bn_set_2_power(p,power,44);
  // cout << p<< endl;
  // cout << FP_PRIME << endl << RLC_DIG << endl << RLC_FP_BITS << endl << WSIZE << endl << RLC_FP_DIGS << endl << p->used << endl;
  // fp_prime_set_dense(p);



  fp_param_set(NIST_256);
  cout << FP_PRIME << endl << RLC_DIG << endl << RLC_FP_BITS << endl << WSIZE << endl << RLC_FP_DIGS << endl;

  const dig_t * p = fp_prime_get();
  const dig_t * q = fp_prime_get_rdc();
  cout << *p << endl << *q << endl;

  fp_param_print();

  bn_t x;
  bn_null(x);
  bn_new(x);
  fp_prime_get_par(x);
  cout << x->used << endl;

  // cout << "1\n";
  // bn_init(a,10);
  // cout << "2\n";
  // bn_set_bit(a,1,11);
  // cout << "3\n";
  // fp_prime_init();
  // cout << "4\n";
  // fp_prime_set_dense(a);
  // // // clean up
  // cout << "5\n";
  // fp_prime_clean();
  // cout << "6\n";
  // bn_free(p);
  // cout << "7\n";
  core_clean();
  return EXIT_SUCCESS;
}
