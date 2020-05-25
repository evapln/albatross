/* implementation of pvss protocol */
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include <unistd.h>

#include <NTL/matrix.h>
#include <NTL/ZZ_pX.h>


#include "pvss.h"

// structures and variables
bool dist = false;
bool rec = false;

// need to be changed in ZZ, and ZZ_p -> instead of arrays : use vectors
// struct pl_t {
//   int n;
//   int t;
//   int q;
//   int l;
//   int r;
//   double h;
//   double *pk;
//   double *pi;
//   double *sig;
//   int *LDEI;
//   int *sigg_ind;
//   double *sigg;
//   int *DLEQ;
//   double *S;
// };

struct pl_t {
  int n; // number of participants
  int t; // threshold
  int l; // number of secrets
  int r; // number of participants that wanting to reconstruct the secrets
  ZZ q; // order of the group G_q
  ZZ_p h; // generator of G_q
  Vec<ZZ_p> pk; // publi keys
  Vec<ZZ_p> sig; // shamir's shares
  Vec<ZZ_p> sighat; // encrypted shares
  int *LDEI;
  Mat<ZZ_p> sigtilde; // decrypted shares and their index
  int *DLEQ;
  Vec<ZZ_p> S; // secrets reconstructed
};

void prng_init(const unsigned int seed) {
  static bool seed_init = false;
  if (! seed_init) {
    srand(seed);
    seed_init = true;
  }
}

// pl_t *pl_alloc(const int n) {
//   pl_t *pl = (pl_t*)malloc(sizeof(pl_t));
//   if (!pl)
//     return NULL;
//   pl->pk = (double*)malloc(n * sizeof(double));
//   pl->pi = (double*)malloc(n * sizeof(double));
//   pl->sig = (double*)malloc(n * sizeof(double));
//   pl->LDEI = (int*)malloc(n * sizeof(int));
//   pl->sigg_ind = (int*)malloc(n * sizeof(int));
//   pl->sigg = (double*)malloc(n * sizeof(double));
//   pl->DLEQ = (int*)malloc(n * sizeof(int));
//   pl->S = (double*)malloc(n * sizeof(double));
//   if (!pl->pk || !pl->pi || !pl->sig || !pl->LDEI || !pl->sigg_ind || !pl->sigg || !pl->DLEQ || !pl->S) {
//     free(pl->S);
//     free(pl->DLEQ);
//     free(pl->sigg_ind);
//     free(pl->sigg);
//     free(pl->LDEI);
//     free(pl->pk);
//     free(pl->pi);
//     free(pl->sig);
//     free(pl);
//     return NULL;
//   }
//   pl->n = n;
//   pl->t = 0;
//   pl->q = 0;
//   pl->l = 0;
//   pl->r = 0;
//   return pl;
// }

pl_t *pl_alloc(const int n) {
  pl_t *pl = (pl_t*)malloc(sizeof(pl_t));
  if (!pl)
    return NULL;
  pl->LDEI = (int*)malloc(n * sizeof(int));
  pl->DLEQ = (int*)malloc(n * sizeof(int));
  if (!pl->LDEI || !pl->DLEQ) {
    free(pl->DLEQ);
    free(pl->LDEI);
    free(pl);
    return NULL;
  }
  pl->n = n;
  pl->t = 0;
  pl->q = ZZ(0);
  pl->l = 0;
  pl->r = 0;
  pl->pk.SetLength(n);
  pl->sig.SetLength(n);
  pl->sighat.SetLength(n);
  return pl;
}

// void pl_free(pl_t *pl) {
//   if (pl) {
//     free(pl->S);
//     free(pl->DLEQ);
//     free(pl->sigg_ind);
//     free(pl->sigg);
//     free(pl->LDEI);
//     free(pl->pk);
//     free(pl->pi);
//     free(pl->sig);
//     free(pl);
//   }
// }

void pl_free(pl_t *pl) {
  if (pl) {
    free(pl->DLEQ);
    free(pl->LDEI);
    free(pl);
  }
}

// void pl_print(pl_t *pl) {
//   puts("\n\n_________________________ Public Ledger _____________________________");
//   printf("  %d participants\n    q = %d\n", pl->n, pl->q);
//   if (pl->pk) {
//     printf("    pk : ");
//     for (int i = 0; i < pl->n; i++)
//       printf("%.0f ", pl->pk[i]);
//   }
//   if (dist) {
//     printf("\n  %d threshold",pl->t);
//     if (pl->pi) {
//       printf("\n    pi : ");
//       for (int i = 0; i < pl->n; i++)
//         printf("%.0f ", pl->pi[i]);
//     }
//     if (pl->sig) {
//       printf("\n    sig : ");
//       for (int i = 0; i < pl->n; i++)
//         printf("%.0f ", pl->sig[i]);
//     }
//     if (pl->LDEI) {
//       printf("\n      LDEI : ");
//       for (int i = 0; i < pl->n; i++)
//         printf("%d ", pl->LDEI[i]);
//     }
//   }
//   if (rec) {
//     printf("\n  %d parties want to reconstruct",pl->r);
//     if (pl->sigg && pl->sigg[0] && pl->sigg[1]) {
//       printf("\n    sigg : ");
//       for (int i = 0; i < pl->r; i++)
//         printf("[%d:%.0f] ", pl->sigg_ind[i], pl->sigg[i]);
//     }
//     if (pl->DLEQ) {
//       printf("\n      DLEQ : ");
//       for (int i = 0; i < pl->r; i++)
//         printf("%d ", pl->DLEQ[i]);
//     }
//     printf("\n  The %d secrets are : ",pl->l);
//     if (pl->S) {
//       for (int i = 0; i < pl->l; i++)
//         printf("%.0f ", pl->S[i]);
//     }
//   }
//   puts("\n_____________________________________________________________________\n");
// }

void pl_print(pl_t *pl) {
  puts("\n\n___________________________ Public Ledger ___________________________");
  printf("%d participants\n",pl->n);
  cout << "q = " << pl->q << endl;
  cout << "pk : " << pl->pk << endl;
  if (dist) {
    printf("\n%d threshold\n",pl->t);
    cout << "shamir's shares :\n" << pl->sig << endl;
    cout << "encrypted shares :\n" << pl->sighat << endl;
    if (pl->LDEI) {
      printf("      LDEI : ");
      for (int i = 0; i < pl->n; i++)
        printf("%d ", pl->LDEI[i]);
    }
    puts("");
  }
  if (rec) {
    printf("\n%d parties want to reconstruct\n",pl->r);
    cout << "decrypted shares :\n";
    for (int i = 1; i <= pl->r; i++)
      cout << pl->sigtilde(i) << " ";
    if (pl->DLEQ) {
      printf("\n      DLEQ : ");
      for (int i = 0; i < pl->r; i++)
        printf("%d ", pl->DLEQ[i]);
    }
    cout << "\n\nThe " << pl->l << " secrets are :\n" << pl->S << endl;
  }
  puts("_____________________________________________________________________\n");
}

void generator(ZZ_p& g, const ZZ& q) {
  // cout << "generator : firt row -> p\t" << ZZ_p::modulus() << endl;
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

// float *setup(const int q, const double h, pl_t *pl) {
//   if (!pl)
//     return NULL;
//   prng_init(time(NULL) + getpid());
//   float *sk = (float*)malloc((pl->n)*sizeof(int));
//   if (!sk)
//     return NULL;
//   float s;
//   for (int i = 0; i < pl->n; i++) {
//     s = rand() % q;
//     while (s == 0)
//       s = rand() % q;
//     sk[i] = s;
//     pl->pk[i] = pow(h, sk[i]);
//     // pow(pl->pk[i], h, sk[i]);
//   }
//   pl->q = q;
//   pl->h = h;
//   return sk;
// }

pl_t *setup(Vec<ZZ_p>& sk, const int n, const ZZ& q, const ZZ& p, const ZZ_p& h) {
  // cout << "setup : firts row -> q\t" << ZZ_p::modulus() << endl;
  ZZ_p s;
  sk.SetLength(n);
  for (int i = 0; i < n; i++) {
    random(s);
    while (IsZero(s))
      random(s);
    sk[i] = s;
  }
  ZZ r;
  ZZ_pPush push(p);
  // cout << "setup : after push -> p\t" << ZZ_p::modulus() << endl;
  pl_t * pl = pl_alloc(n);
    if (!pl)
      return NULL;
  for (int i = 0; i < n; i++) {
    r = rep(sk[i]);
    power(pl->pk[i], h, r);
  }
  pl->n = n;
  pl->q = q;
  pl->h = h;
  return pl;
}

// void distribution(const int l, const int t, pl_t *pl) {
//   if (!pl || t < 1 || t > pl->n)
//     return;
//   int deg = t + l;
//   int *p = (int*)rand_poly(deg, pl->q);
//   if (p) {
//     pl->l = l;
//     pl->t = t;
//     printf("poly = %d", p[0]);
//     for (int i = 1; i < deg; i++) {
//       printf(" + %d*x^%d",p[i],i);
//     }
//     puts("");
// // print secrets for verification //////////////////////////
//     puts("\nSecrets :");
//     for (int i = 0; i < l; i++) {
//       double ap = apply_poly(p,deg,-i, pl->q);
//       printf("S%d = h^%f = %f\n", i, ap, pow(pl->h,ap));
//     }
// ////////////////////////////////////////////////////////////
//     for (int i = 0; i < pl->n; i++) {
//       double ap = apply_poly(p, deg, i+1, pl->q);
//       pl->pi[i] = ap;
//       pl->sig[i] = pow(pl->pk[i],ap);
//       pl->LDEI[i] = 0; // PREUVE A FAIRE
//     }
//     dist = true;
//     free(p);
//   }
// }

void distribution(const int l, const int t, pl_t *pl) {
  // cout << "distribution : first row -> q\t" << ZZ_p::modulus() << endl;
  if (!pl || t < 1 || t > pl->n)
    return;
  int deg = t + l;
  ZZ_pX P;
  random(P,deg);
  pl->l = l;
  pl->t = t;
  Vec<ZZ_p> s;
  s.SetLength(pl->n+l);
  ZZ_p tmp;
  ZZ repzz;
  // computation of s0,...,s(l-1) (secrets) and sig1,...,sign (shamir's shares)
  for (int i = -l+1; i <= pl->n; i++) {
    tmp = ZZ_p(i);
    eval(s[i+l-1], P, tmp);
  }
  ZZ p = 2 * pl->q + 1;
  ZZ_pPush push(p);
  // cout << "distribution : after push -> p\t" << ZZ_p::modulus() << endl;
  // attribution of the shamir's shares their values and computation of encrypted shares
  for (int i = 0; i < pl->n; i++) {
    pl->sig[i] = s[i+l];
    repzz = rep(s[i+l]);
    power(pl->sighat[i],pl->pk[i],repzz);
    pl->LDEI[i] = 0; // PREUVE A FAIRE
  }
  // print secrets for verification //////////////////////////
  puts("\nSecrets :");
  // computation of secrets
  for (int i = l-1; i >= 0; i--) {
    repzz = rep(s[i]);
    power(tmp,pl->h,repzz);
    cout << "S" << i << " = h^" << repzz << " = " << tmp << endl;
  }
  dist = true;
  // cout << "distribution : end -> p\t" << ZZ_p::modulus() << endl;
}

// int lambda(const int i, const int j, const int *tab, const int len) {
//   double num, den;
//   num = 1;
//   den = 1;
//   for (int k = 0; k < len; k++) {
//     double m = tab[k];
//     if (m != i) {
//       num *= (-j-m);
//       den *= (i-m);
//     }
//   }
//   return num/den;
// }

// void lambda(Mat<ZZ_p>& lambs, int t, pl_t *pl) {
//   lambs.SetDims(t,pl->l);
//   ZZ num;
//   ZZ den;
//   ZZ di;
//   ZZ_p tmp;
//   for (int j = 0; j < pl->l; j++) {
//     for (int i = 1; i <= t; i++) {
//       num = ZZ(1); den = ZZ(1);
//       for (int m = 1; m <= t; m++)
//         if (m != i) {
//           sub(tmp,- ZZ_p(j),pl->sigtilde(m,1));
//           mul(num, num, rep(tmp));
//           sub(tmp, pl->sigtilde(i,1), pl->sigtilde(m,1));
//           mul(den, den, rep(tmp));
//         }
//       div(di,num,den);
//       conv(lambs(i,j+1),di);
//     }
//   }
// }

void lambda(Mat<ZZ_p>& lambs, int t, pl_t *pl) {
  lambs.SetDims(t,pl->l);
  ZZ_p num;
  ZZ_p den;
  ZZ_p mu;
  ZZ_p tmp;
  ZZ_p invden;
  for (int j = 0; j < pl->l; j++) {
    for (int i = 1; i <= t; i++) {
      num = ZZ_p(1); den = ZZ_p(1);
      for (int m = 1; m <= t; m++)
        if (m != i) {
          sub(tmp,- ZZ_p(j),pl->sigtilde(m,1));
          mul(num, num, tmp);
          sub(tmp, pl->sigtilde(i,1), pl->sigtilde(m,1));
          mul(den, den, tmp);
        }
      inv(invden,den);
      mul(mu,num,invden);
      conv(lambs(i,j+1),mu);
    }
  }
}

// void reconstruction(Vec<ZZ>& S, const int r, pl_t *pl) {
//   if (!pl) // error : the public ledger doesn't exist
//     return;
//   int t = pl->n - pl->t;
//   if (r < t) // error : not enough parts
//     return;
//   S.SetLength(pl->l);
//   // verification of proof DLEQ
//   double num, den;
//   float s;
//   for (int j = 0; j < pl->l; j++) {
//     s = 1;
//     for (int i = 0; i < t; i++) {
//       num = 1; den = 1;
//       for (int m = 0; m < t; m++)
//         if (m != i) {
//           num *= - j - pl->sigg_ind[m];
//           den *= (pl->sigg_ind[i] - pl->sigg_ind[m]);
//         }
//       // the following line doesn't work because pow can't take a ZZ in power
//       s *= pow(pl->sigg[i],num/den);
//     }
//     S[j] = s;
//     pl->S[j] = 0; // need to declare pl->S another way
//   }
//   pl->r = r;
//   rec = true;
// }

void reconstruction(const int r, pl_t *pl) {
  // cout << "reconstruction : first row -> p\t" << ZZ_p::modulus() << endl;
  if (!pl) // error : the public ledger doesn't exist
    return;
  int t = pl->n - pl->t;
  if (r < t) // error : not enough parts
    return;
  // verification of proof DLEQ
  pl->S.SetLength(pl->l);
  ZZ_pPush push(pl->q);
  // cout << "reconstruction : after push 1 -> q\t" << ZZ_p::modulus() << endl;
  Mat<ZZ_p> lambs;
  lambda(lambs, t, pl);
  // cout << "lambdas: " << lambs << endl;
  ZZ lamb;
  ZZ_p::init(2 * pl->q + 1);
  // cout << "reconstruction : after push 2 -> p\t" << ZZ_p::modulus() << endl;
  ZZ_p tmp;
  for (int j = 0; j < pl->l; j++) {
    pl->S[j] = ZZ_p(1);
    for (int i = 1; i <= t; i++) {
      lamb = rep(lambs(i,j+1));
      power(tmp,pl->sigtilde(i,2),lamb);
      mul(pl->S[j], pl->S[j], tmp);
    }
  }
  pl->r = r;
  rec = true;
}

int main(void) {
  int n = 20;
  ZZ q;
  // q = ZZ(7);
  GenGermainPrime(q,10);
  ZZ p = 2 * q + 1;
  // cout << "p = " << p << "\tq = " << q << endl;
  int t = 6;
  int l = n-2*t;
  ZZ_p::init(p);
  // cout << "main : after init 1 -> p\t" << ZZ_p::modulus() << endl;
  ZZ_p g;
  generator(g,p);
  ZZ_p h;
  power(h,g,2);
  // prng_init(time(NULL) + getpid());
  // double h = 3; //generator of F_7
  ZZ_p::init(q);
  // cout << "main : after inti 2 -> q\t" << ZZ_p::modulus() << endl;
  Vec<ZZ_p> sk;
  pl_t *pl = setup(sk,n,q,p,h);
  // cout << "main : after setup -> q\t" << ZZ_p::modulus() << endl;
  // cout << "scrape" << endl;
  if (!pl) {
    puts("oups");
    return EXIT_FAILURE;
  }
  cout << endl << "sk : " << sk << endl;

  // pl_print(pl);
  distribution(l,t,pl);
  // cout << "main : after distribution -> q\t" << ZZ_p::modulus() << endl;
  // pl_print(pl);

  // choice of the r participants who want to recover the secret vector
  int tab[n];
  int len = n;
  int r = n - t;
  int choice[r];
  pl->sigtilde.SetDims(r,2);
  Vec<ZZ_p> invsk;
  invsk.SetLength(r);
  ZZ_p tmp;
  prng_init(time(NULL) + getpid());
  for (int i = 0; i < len; i++)
    tab[i] = i;
  for (int i = 1; i <= r; i++) {
    int ind = rand() % len;
    int in = tab[ind];
    pl->sigtilde(i,1) = ZZ_p(in + 1);
    inv(invsk[i],sk[in]);
    choice[i] = in;
    len--;
    for (int j = ind; j < len; j++)
      tab[j] = tab[j+1];
  }
  ZZ_p::init(p);
  // cout << "distribution : after init 3 -> p\t" << ZZ_p::modulus() << endl;
  for (int i = 1; i <= r; i++) {
    power(tmp,pl->sighat[choice[i]],rep(invsk[i]));
    pl->sigtilde(i,2) = tmp;
    pl->DLEQ[i] = 0;
  }
  // cout << "choice done" << endl;

  // cout << "main : after choice -> p\t" << ZZ_p::modulus() << endl;

  // pl_print(pl);
  reconstruction(r, pl);
  // cout << "main : after reconstruction -> p\t" << ZZ_p::modulus() << endl;
  pl_print(pl);

// les decrypted shares sont correctes
  // ZZ_p verif;
  // ZZ re;
  // cout << "real decrypted shares" << endl;
  // for (int i = 0; i < pl->n; i++) {
  //   re = rep(pl->sig[i]);
  //   power(verif,pl->h,re);
  //   cout << "[" << i+1 << " " << verif << "] ";
  // }

  pl_free(pl);
  printf("\n\n---------------------------------------------------------------------\n");
  // cout << "main : end -> p\t" << ZZ_p::modulus() << endl;


  return EXIT_SUCCESS;
}


/* TODO :
- implémentation figure 6 dans scrape++
- nettoyer le code */
