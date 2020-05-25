// some tries of shamir secret sharing

// compile with :
// g++ -c -Wall -Wextra  draft.c
// g++ -o draft draft.o -lntl -lgmp -pthread



#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <unistd.h>

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/ZZ_pX.h>

using namespace std;
using namespace NTL;

void prng_init(const unsigned int seed) {
  static bool seed_init = false;
  if (! seed_init) {
    srand(seed);
    seed_init = true;
  }
}

int mod(const int n, const int q) {
  int r = n % q;
  return r < 0 ? r + q :r;
}

double apply_poly(const int *poly, const int deg, const double val, const int q) {
  double res = 0;
  long mul = 1;
  for (int i = 0; i <= deg; i++) {
    // res += poly[i] * mul;
    res = mod(res + mod(poly[i] * mul,q),q);
    // mul *= val;
    mul = mod(mul * val, q);
  }
  return res;
}

int *rand_poly(const int deg, const int q) {
  if (deg < 0) // non valid degree
    return NULL;
  prng_init(time(NULL) + getpid());
  int *p = (int*)malloc((deg+1)*sizeof(int));
  if (!p) // error of allocation
    return NULL;
  for (int i = 0; i <= deg; i++)
    p[i] = rand() % q;
  return p;
}

// sharing of one secret : S = p(0)
double *sharing(const int t, const int n, const int q) {
  int *p = rand_poly(t-1,q);
  if (!p)
    return NULL;
  double secret = apply_poly(p,t-1,0,q);
  printf("poly = %d", p[0]);
  for (int i = 1; i < t; i++) {
    printf(" + %dx^%d",p[i],i);
  }
  printf("\nsecret : %.0f\n", secret);
  double *shares = (double*)malloc(n * sizeof(double));
  if (!shares)
    return NULL;
  for (int i = 0; i < n; i++)
    shares[i] = apply_poly(p,t-1,i+1,q);
  free(p);
  return shares;
}

// sharing of l secrets : Si = p(-i) for i in {0,l-1}
double *sharing_multi(const int t, const int n, const int l, const int q) {
  int *p = rand_poly(t-1,q);
  double *secret = (double*)malloc(l*sizeof(double));
  if (!secret)
    return NULL;
  for (int i = 0; i < l; i ++)
    secret[i] = apply_poly(p,t-1,-i,q);
  printf("\npoly = %d", p[0]);
  for (int i = 1; i < t; i++)
    printf(" + %d*x^%d",p[i],i);
  printf("\nsecrets: ");
  for (int i = 0; i < l; i++)
    printf("S%d = %.0f\t", i, secret[i]);
  double *shares = (double*)malloc(n * sizeof(double));
  for (int i = 0; i < n; i++)
    shares[i] = apply_poly(p,t-1,i+1,q);
  free(p);
  free(secret);
  return shares;
}

// reconstruction of one secret
double reconstruction(double **shares, const int nb_shares, const int t,
                   const int q) {
  if (nb_shares < t || !shares)
    return -1;
  double s = 0, num, den;
  for (int j = 0; j < nb_shares; j++) {
    num = 1; den = 1;
    for (int m = 0; m < nb_shares; m++) {
      if (m != j) {
        num *= -shares[m][0];
        den *= (shares[j][0] - shares[m][0]);
      }
    }
    s = mod(s+shares[j][1] * num / den,q);
    // printf("\nyj = %.0f, num = %.0f, den = %.0f, sum = %.0f",shares[j][1] ,num ,den, shares[j][1] * num / den);
  }
  return s;
}

// reconstruction of l secrets
// (problems when two many secrets -> need a greater size of variables)
double *reconstruction_multi(double **shares, const int nb_shares, const int l, const int t,
                   const int q) {
  if (nb_shares < t)
    return NULL;
  double *s = (double*)malloc(l*sizeof(double));
  for (int j = 0; j < l; j++)
    s[j] = 0;
  double num, den;
  for (int j = 0; j < l; j++) {
    for (int i = 0; i < t; i++) {
      num = 1; den = 1;
      for (int m = 0; m < t; m++) {
        if (m != i) {
          num *= -j-shares[m][0];
          den *= (shares[i][0] - shares[m][0]);
        }
      }
      s[j] = mod(s[j]+mod(shares[i][1] * mod(num / den,q),q),q);
    }
  }
  return s;
}

// reconstruction of l secrets with the library NTL
void reconstruction_multiZZ(Vec<ZZ>& s, double **shares, const int nb_shares, const int l, const int t,
                   const int q) {
  if (nb_shares < t)
    return;
  s.SetLength(l);
  for (int j = 0; j < l; j++)
    s[j] = 0;
  ZZ num;
  ZZ den;
  for (int j = 0; j < l; j++) {
    for (int i = 0; i < t; i++) {
      num = ZZ(1); den = ZZ(1);
      for (int m = 0; m < t; m++)
        if (m != i) {
          num *= -j-shares[m][0];
          den *= (shares[i][0] - shares[m][0]);
        }
      s[j] = (s[j] + (shares[i][1] * (num / den) % q ) % q) % q;
    }
  }
}



// main for reconstruction_multiZZ
// int main(void) {
//   int n = 60;
//   int t = 16;
//   int q = 151;
//   int r = 16;
//   int l = n-2*t;
//   double *share = sharing_multi(t, n, l, q);
//   double **shares = (double**)malloc(r*sizeof(double*));
//   for (int i = 0; i < r; i++)
//     shares[i] = (double*)malloc(2*sizeof(double));
//   for (int i = 0; i < r; i++) {
//     shares[i][0] = i+1;
//     shares[i][1] = share[(int)shares[i][0]-1];
//   }
//   for (int i = 0; i < r; i++)
//     printf("[%.0f,%.0f] ",shares[i][0],shares[i][1]);
//   Vec<ZZ> s;
//   reconstruction_multiZZ(s, shares, r, l, t, q);
//   cout << std::endl << s << std::endl;
//   puts("");
//   free(share);
//   for (int i = 0; i < r; i++)
//     free(shares[i]);
//   free(shares);
// }

// // main for reconstruction_multi
// int main(void) {
//   int n = 60;
//   int t = 16;
//   int q = 151;
//   int r = 16;
//   int l = n-2*t;
//   double *share = sharing_multi(t, n, l, q);
//   double **shares = (double**)malloc(r*sizeof(double*));
//   for (int i = 0; i < r; i++)
//     shares[i] = (double*)malloc(2*sizeof(double));
//   for (int i = 0; i < r; i++) {
//     shares[i][0] = i+1;
//     shares[i][1] = share[(int)shares[i][0]-1];
//   }
//   for (int i = 0; i < r; i++)
//     printf("[%.0f,%.0f] ",shares[i][0],shares[i][1]);
//   double *s = reconstruction_multi(shares, r, l, t, q);
//   for (int i = 0; i < l; i++)
//     printf("\nS%d = %.0f",i,s[i]);
//   puts("");
//   free(share);
//   for (int i = 0; i < r; i++)
//     free(shares[i]);
//   free(shares);
//   free(s);
// }

// // main for reconstruct
// int main(void) {
//   int n = 11;
//   int t = 6;
//   int q = 17;
//   int r = 7;
//   double *share = sharing(t, n, q);
//   if(!share) {printf("error\n"); return EXIT_FAILURE;}
//   double **shares = malloc(r*sizeof(double*));
//   if(!shares) {printf("error2\n"); return EXIT_FAILURE;}
//   for (int i = 0; i < r; i++)
//     shares[i] = malloc(2*sizeof(double));
//   for (int i = 0; i < r; i++) {
//     shares[i][0] = i+1;
//     shares[i][1] = share[(int)shares[i][0]-1];
//   }
//   for (int i = 0; i < r; i++)
//     printf("[%.0f,%.0f] ",shares[i][0],shares[i][1]);
//   double s = reconstruction(shares, r, t, q);
//   printf("\nsecret = %.0f\n",s);
//   free(share);
//   for (int i = 0; i < r; i++)
//     free(shares[i]);
//   free(shares);
// }

////////////////////////////////////////////////////////////////////////////////
//////////////////////// Primes and all of that ////////////////////////////////

void generator(ZZ_p& g, const ZZ& q) {
  ZZ p;
  ZZ_p po;
  for (int i = 2; i < 2*q+1; i++) {
    power(po,ZZ_p(i),2);
    if (po == 1)
      continue;
    power(po,ZZ_p(i),q);
    if (IsOne(po))
      continue;
    g = ZZ_p(i);
    return;
  }
  return;
}

void change_modulus(ZZ& p) {
  // ZZ_p::init(q);
  cout << ZZ_p::modulus() << endl;
  // ZZ_p::init(p);
  ZZ_pPush push(p);
  cout << ZZ_p::modulus() << endl;
}
// test sophie gremain prime with NTL
int main(void) {
  // ZZ q;
  // GenGermainPrime(q,10);
  // cout << "q = " << q << endl;
  // ZZ_p::init(2*q+1);
  // // ZZ_p g;
  // // generator(g,q);
  // // cout << "g = " << g << endl;
  // Mat<ZZ_p> A;
  // int n = 3, m = 4;
  // A.SetDims(n,m);
  // cout << "ok\n";
  // for (int i = 1; i <= n; i++)
  //   for (int j = 1; j <= m; j++)
  //     A(i,j) = m*(i-1)+j;
  // cout << "ok\n";
  // cout << A << endl;
  // ZZ_p::init(ZZ(11));
  // ZZ_pX p;
  // random(p,4);
  // cout << p << endl;
  // ZZ_p a;
  // ZZ_p b;
  // cin >> a;
  // eval(b,p,a);
  // cout << b << endl;
  //
  // ZZ_p in;
  // for (int i = 1; i < 11; i++) {
  //   inv(in,ZZ_p(i));
  //   cout << i << " -> " << in << endl;
  // }

  ZZ p;
  ZZ q;
  // ZZ mod;
  GenGermainPrime(q,10);
  // q  = 7;
  p = 2*q+1;
  // q = 2;
  cout << "q = " << q << "\np = " << p << endl;
  ZZ_p::init(q);
  cout << "mod in main -> " << ZZ_p::modulus() << endl;
  change_modulus(p);
  cout << "mod in main -> " << ZZ_p::modulus() << endl;
  Vec<ZZ_p> v;
  v.SetLength(5);
  for (int i = 0; i < 5; i++)
    v[i] = i;
  cout << v << endl;
  return EXIT_SUCCESS;
}
