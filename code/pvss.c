/* implementation of pvss protocol */
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include <unistd.h>


#include "pvss.h"

// structures and variables
bool dist = false;
bool rec = false;

// need to be changed in ZZ, and ZZ_p -> instead of arrays : use vectors
struct pl_t {
  int n;
  int t;
  int q;
  int l;
  int r;
  double h;
  double *pk;
  double *pi;
  double *sig;
  int *LDEI;
  int *sigg_ind;
  double *sigg;
  int *DLEQ;
  double *S;
};

pl_t *pl_alloc(const int n) {
  pl_t *pl = (pl_t*)malloc(sizeof(pl_t));
  if (!pl)
    return NULL;
  pl->pk = (double*)malloc(n * sizeof(double));
  pl->pi = (double*)malloc(n * sizeof(double));
  pl->sig = (double*)malloc(n * sizeof(double));
  pl->LDEI = (int*)malloc(n * sizeof(int));
  pl->sigg_ind = (int*)malloc(n * sizeof(int));
  pl->sigg = (double*)malloc(n * sizeof(double));
  pl->DLEQ = (int*)malloc(n * sizeof(int));
  pl->S = (double*)malloc(n * sizeof(double));
  if (!pl->pk || !pl->pi || !pl->sig || !pl->LDEI || !pl->sigg_ind || !pl->sigg || !pl->DLEQ || !pl->S) {
    free(pl->S);
    free(pl->DLEQ);
    free(pl->sigg_ind);
    free(pl->sigg);
    free(pl->LDEI);
    free(pl->pk);
    free(pl->pi);
    free(pl->sig);
    free(pl);
    return NULL;
  }
  pl->n = n;
  pl->t = 0;
  pl->q = 0;
  pl->l = 0;
  pl->r = 0;
  return pl;
}

void pl_free(pl_t *pl) {
  if (pl) {
    free(pl->S);
    free(pl->DLEQ);
    free(pl->sigg_ind);
    free(pl->sigg);
    free(pl->LDEI);
    free(pl->pk);
    free(pl->pi);
    free(pl->sig);
    free(pl);
  }
}

void pl_print(pl_t *pl) {
  puts("\n\n_________________________ Public Ledger _____________________________");
  printf("  %d participants\n    q = %d\n", pl->n, pl->q);
  if (pl->pk) {
    printf("    pk : ");
    for (int i = 0; i < pl->n; i++)
      printf("%.0f ", pl->pk[i]);
  }
  if (dist) {
    printf("\n  %d threshold",pl->t);
    if (pl->pi) {
      printf("\n    pi : ");
      for (int i = 0; i < pl->n; i++)
        printf("%.0f ", pl->pi[i]);
    }
    if (pl->sig) {
      printf("\n    sig : ");
      for (int i = 0; i < pl->n; i++)
        printf("%.0f ", pl->sig[i]);
    }
    if (pl->LDEI) {
      printf("\n      LDEI : ");
      for (int i = 0; i < pl->n; i++)
        printf("%d ", pl->LDEI[i]);
    }
  }
  if (rec) {
    printf("\n  %d parties want to reconstruct",pl->r);
    if (pl->sigg && pl->sigg[0] && pl->sigg[1]) {
      printf("\n    sigg : ");
      for (int i = 0; i < pl->r; i++)
        printf("[%d:%.0f] ", pl->sigg_ind[i], pl->sigg[i]);
    }
    if (pl->DLEQ) {
      printf("\n      DLEQ : ");
      for (int i = 0; i < pl->r; i++)
        printf("%d ", pl->DLEQ[i]);
    }
    printf("\n  The %d secrets are : ",pl->l);
    if (pl->S) {
      for (int i = 0; i < pl->l; i++)
        printf("%.0f ", pl->S[i]);
    }
  }
  puts("\n_____________________________________________________________________\n");
}

int mod(const int n, const int q) {
  int r = n % q;
  return r < 0 ? r + q :r;
}


void prng_init(const unsigned int seed) {
  static bool seed_init = false;
  if (! seed_init) {
    srand(seed);
    seed_init = true;
  }
}

double apply_poly(const int *poly, const int deg, const int val, const int q) {
  long mul = 1;
  double res = 0;
  for (int i = 0; i <= deg; i++) {
    res = mod(res + mod(poly[i] * mul,q),q);
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

float *setup(const int q, const double h, pl_t *pl) {
  if (!pl)
    return NULL;
  prng_init(time(NULL) + getpid());
  float *sk = (float*)malloc((pl->n)*sizeof(int));
  if (!sk)
    return NULL;
  float s;
  for (int i = 0; i < pl->n; i++) {
    s = rand() % q;
    while (s == 0)
      s = rand() % q;
    sk[i] = s;
    pl->pk[i] = pow(h, sk[i]);
    // pow(pl->pk[i], h, sk[i]);
  }
  pl->q = q;
  pl->h = h;
  return sk;
}

void distribution(const int l, const int t, pl_t *pl) {
  if (!pl || t < 1 || t > pl->n)
    return;
  int deg = t + l;
  int *p = (int*)rand_poly(deg, pl->q);
  if (p) {
    pl->l = l;
    pl->t = t;
    printf("poly = %d", p[0]);
    for (int i = 1; i < deg; i++) {
      printf(" + %d*x^%d",p[i],i);
    }
    puts("");
// print secrets for verification //////////////////////////
    puts("\nSecrets :");
    for (int i = 0; i < l; i++) {
      double ap = apply_poly(p,deg,-i, pl->q);
      printf("S%d = h^%f = %f\n", i, ap, pow(pl->h,ap));
    }
////////////////////////////////////////////////////////////
    for (int i = 0; i < pl->n; i++) {
      double ap = apply_poly(p, deg, i+1, pl->q);
      pl->pi[i] = ap;
      pl->sig[i] = pow(pl->pk[i],ap);
      pl->LDEI[i] = 0; // PREUVE A FAIRE
    }
    dist = true;
    free(p);
  }
}

int lambda(const int i, const int j, const int *tab, const int len) {
  double num, den;
  num = 1;
  den = 1;
  for (int k = 0; k < len; k++) {
    double m = tab[k];
    if (m != i) {
      num *= (-j-m);
      den *= (i-m);
    }
  }
  return num/den;
}

 void reconstruction(Vec<ZZ>& S, const int r, pl_t *pl) {
  if (!pl) // error : the public ledger doesn't exist
    return;
  int t = pl->n - pl->t;
  if (r < t) // error : not enough parts
    return;
  S.SetLength(pl->l);
  // verification of proof DLEQ
  double num, den;
  float s;
  for (int j = 0; j < pl->l; j++) {
    s = 1;
    // cout <<  j << " : " << s << endl;
    for (int i = 0; i < t; i++) {
      num = 1; den = 1;
      for (int m = 0; m < t; m++)
        if (m != i) {
          num *= - j - pl->sigg_ind[m];
          den *= (pl->sigg_ind[i] - pl->sigg_ind[m]);
        }
      //// cout << "\nnum : " << num << "\nden : "<< den << "\ndiv : " << num/den << endl;
      //// int lamb = lambda(pl->sigg_ind[i],-j,pl->sigg_ind,r);
      //// sum += (s[j] + (pl->pi[i] * (num / den) % q ) % q) % q;
      //// operator*=(s,pow(pl->sigg[i],lambda(pl->sigg_ind[i],j,pl->sigg_ind,r)));
      //// printf("sigg = %f, j = %d, i = %d\n",pl->sigg[i], j, pl->sigg_ind[i]);

      // the following line doesn't work because pow can't take a ZZ in power
      s *= pow(pl->sigg[i],num/den);
      // cout << j << " : " << num/den << " " << pow(pl->sigg[i],num/den) << " " << s << endl;
    }
    S[j] = s;
    // printf("p(%d) = %f\t",-j,sum);
    pl->S[j] = 0;
  }
  pl->r = r;
  rec = true;
  // return NULL;
}


int main(void) {
  int n = 6;
  int q = 7;
  int t = 2;
  int l = n-2*t;
  prng_init(time(NULL) + getpid());
  double h = 3; //generator of F_7
  pl_t *pl = pl_alloc(n);
  if (!pl) {
    puts("oups");
    return EXIT_FAILURE;
  }
  float *sk = setup(q,h,pl);
  if (!sk) {
    puts("error");
    pl_free(pl);
    return EXIT_FAILURE;
  }
  printf("sk : ");
  for (int i = 0; i < n; i++)
    printf("%f ",sk[i]);
  puts("\n");

  pl_print(pl);
  distribution(l,t,pl);
  pl_print(pl);

  // choice of the r participants who want to recover the secret vector
  int tab[n];
  int len = n;
  int r = n - t;
  for (int i = 0; i < len; i++)
    tab[i] = i;
  for (int i = 0; i < r; i++) {
    int ind = rand() % len;
    int in = tab[ind];
    pl->sigg_ind[i] = in + 1;
    pl->sigg[i] = trunc(pow(pl->sig[in],1/sk[in]));
    pl->DLEQ[i] = 0;
    len--;
    for (int j = ind; j < len; j++)
      tab[j] = tab[j+1];
  }

  printf("\n%d participants want to recover the secret vector: \n", r);
  for (int i = 0; i < r; i++)
    printf(" (%d,%.0f)\n",pl->sigg_ind[i], pl->sigg[i]);

  Vec<ZZ> s;
  reconstruction(s, r, pl);
  pl_print(pl);

  cout << "bad secrets recovered: " << s <<endl;

  free(sk);
  pl_free(pl);
  printf("\n\n---------------------------------------------------------------------\n");


  return EXIT_SUCCESS;
}


/* TODO :
- voir d'où vient l'erreur dans setup !!!!!!!!!
  -> remplacer tableaux par vecteurs
  -> se servir du corps dans NTL
- créer une fonction renvoyant un générateur d'un groupe multiplicatif
- régler problème de reconstruction
  -> voir si sigg correspondent à h^sig
- implémentation figure 6 dans scrape++
- ajouter une fonction alloc qui permet de n'allouer que l'espace nécessaire
- nettoyer le code */
