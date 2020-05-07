/* implementation of pvss protocol */
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include <unistd.h>

#include <NTL/ZZ.h>

#include "pvss.h"

// structures and variables
bool dist = false;
bool rec = false;

struct pl_t {
  int n;
  int t;
  int q;
  int l;
  int r;
  int h;
  int *pk;
  int *sig;
  int *LDEI;
  int **sigg;
  int *DLEQ;
  int *S;
};

pl_t *pl_alloc(int n) {
  pl_t *pl = (pl_t*)malloc(sizeof(pl_t));
  if (!pl)
    return NULL;
  pl->pk = (int*)malloc(n * sizeof(int));
  pl->sig = (int*)malloc(n * sizeof(int));
  pl->LDEI = (int*)malloc(n * sizeof(int));
  pl->sigg = (int**)malloc(2 * sizeof(int*));
  pl->DLEQ = (int*)malloc(n * sizeof(int));
  pl->S = (int*)malloc(n * sizeof(int));
  if (!pl->pk || !pl->sig || !pl->LDEI || !pl->sigg || !pl->DLEQ || !pl->S) {
    free(pl->S);
    free(pl->DLEQ);
    free(pl->sigg);
    free(pl->LDEI);
    free(pl->pk);
    free(pl->sig);
    free(pl);
    return NULL;
  }
  pl->sigg[0] = (int*)malloc(n*sizeof(int));
  pl->sigg[1] = (int*)malloc(n*sizeof(int));
  if (!pl->sigg[0] || !pl->sigg[0]) {
    free(pl->S);
    free(pl->DLEQ);
    free(pl->sigg);
    free(pl->LDEI);
    free(pl->pk);
    free(pl->sig);
    free(pl->sigg[0]);
    free(pl->sigg[1]);
    free(pl->sigg);
    free(pl);
    return NULL;
  }
  pl->n = n;
  pl->t = 0;
  pl->q = 0;
  pl->l = 0;
  pl->r = 0;
  pl->h = 0;
  return pl;
}

void pl_free(pl_t *pl) {
  if (pl) {
    free(pl->S);
    free(pl->DLEQ);
    free(pl->sigg[0]);
    free(pl->sigg[1]);
    free(pl->sigg);
    free(pl->LDEI);
    free(pl->pk);
    free(pl->sig);
    free(pl);
  }
}

void pl_print(pl_t *pl) {
  puts("\nPublic Ledger :");
  printf("%d participants\n%d threshold\nq = %d\n", pl->n, pl->t, pl->q);
  if (pl->pk) {
    printf("pk : ");
    for (int i = 0; i < pl->n; i++)
      printf("%d ", pl->pk[i]);
  }
  if (dist) {
    if (pl->sig) {
      printf("\nsig : ");
      for (int i = 0; i < pl->n; i++)
        printf("%d ", pl->sig[i]);
    }
    if (pl->LDEI) {
      printf("\nLDEI : ");
      for (int i = 0; i < pl->n; i++)
        printf("%d ", pl->LDEI[i]);
    }
  }
  if (rec) {
    printf("\n%d parties want to reconstruct\n",pl->r);
    if (pl->sigg && pl->sigg[0] && pl->sigg[1]) {
      printf("\nsigg : ");
      for (int i = 0; i < pl->r; i++)
        printf("[%d:%d] ", pl->sigg[0][i], pl->sigg[1][i]);
    }
    if (pl->DLEQ) {
      printf("\nDLEQ : ");
      for (int i = 0; i < pl->r; i++)
        printf("%d ", pl->DLEQ[i]);
    }
    printf("\nThe %d secrets are\n",pl->l);
    if (pl->S) {
      for (int i = 0; i < pl->l; i++)
        printf("%d ", pl->S[i]);
    }
  }
  puts("");
}

int mod(int n, int q) {
  int r = n % q;
  return r < 0 ? r + q :r;
}

int inv(int a, int q) {
  for (int b = 1; b < q; b++)
    if (mod(a*b,q) == 1)
      return b;
  return 0;
}

// initialisation of randomness -> call with seed = time(NULL) + getpid()
void prng_init(unsigned int seed) {
  static bool seed_init = false;
  if (! seed_init) {
    srand(seed);
    seed_init = true;
  }
}

int apply_poly(int *poly, int deg, int val, int q) {
  int res = 0;
  long mul = 1;
  for (int i = 0; i < deg; i++) {
    res = mod(res + mod(poly[i] * mul,q),q);
    mul = mod(mul * val, q);
  }
  return res;
}

// random polynomial in Z_q[X] of degree deg
int *rand_poly(int deg, int q) {
  if (deg < 0) // non valid degree
    return NULL;
  prng_init(time(NULL) + getpid());
  int *p = (int*)malloc(deg*sizeof(int));
  if (!p) // error of allocation
    return NULL;
  for (int i = 0; i < deg; i++)
    p[i] = rand() % q;
  return p;
}

// setup return the secretkey
int *setup(int q, int h, pl_t *pl) {
  if (!pl)
    return NULL;
  prng_init(time(NULL) + getpid());
  int *sk = (int*)malloc(pl->n*sizeof(int));
  if (!sk)
    return NULL;
  for (int i = 0; i < pl->n; i++) {
    sk[i] = rand() % q;
    pl->pk[i] = mod(h^sk[i], q);
  }
  pl->q = q;
  pl->h = h;
  return sk;
}


// n participants, t threshold, l = n-2t-1
// return the list of the secrets ( one secret for one participant)
void distribution(int l, int t, pl_t *pl) {
  if (!pl || t < 1 || t > pl->n)
    return;
  int deg = t + l;
  int *p = (int*)rand_poly(deg, pl->q);
  if (p) {
    pl->l = l;
    pl->t = t;
    puts("poly ");
    for (int i = 0; i < deg; i++)
      printf("%d ",p[i]);
    puts("");
    // print secrets for verification
    puts("\nSecrets :");
    for (int i = 0; i < l; i++) {
      int ap = apply_poly(p,deg,-i, pl->q);
      printf("S%d = h^%d = %d\n", i, ap, mod((pl->h)^ap, pl->q));
    }
    for (int i = 0; i < pl->n; i++) {
      int ap = apply_poly(p, deg, i+1, pl->q);
      // printf("%d : %d\n",i+1,ap);
      pl->sig[i] = mod((pl->pk[i])^ap, pl->q);
      pl->LDEI[i] = 0; // PREUVE A FAIRE
    }
    dist = true;
    free(p);
  }
}

int lambda(int i, int j, int *tab, int len) {
  int num = 1, den = 1;
  for (int k = 0; k < len; k++) {
    int m = tab[k];
    if (m != i) {
      num *= (-j-m);
      den *= (i-m);
    }
  }
  return num/den;
}

// reconstruct the secret vector
int *reconstruction(int r, pl_t *pl) {
  if (!pl) // error : the public ledger doesn't exist
    return NULL;
  if (r < pl->n - pl->t) // error : not enough parts
    return NULL;
  // verification of proof DLEQ
  for (int j = 0; j < pl->l; j++) {
    int s = 1;
    for (int i = 0; i < r; i++) {
      s = mod(s * pl->sigg[1][i]^lambda(pl->sigg[0][i],j,pl->sigg[0],r),pl->q);
    }
    pl->S[j] = s;
  }
  pl->r = r;
  rec = true;
  return NULL;
}

// proofs : LDEI DLEQ
////////////////////////////////////////////////////////////////////////////////



/////// main
int main(void) {
  int n = 32;
  int q = 41;
  int t = 10;
  int l = n-2*t;
  prng_init(time(NULL) + getpid());
  int h = 0; //generator of F_q (q prime so every non zero element is a generator)
  while (h == 0)
    h = rand() % q;
  pl_t *pl = (pl_t*)pl_alloc(n);
  if (!pl) {
    puts("oups");
    return EXIT_FAILURE;
  }
  int *sk = (int*)setup(q,h,pl);
  if (!sk) {
    puts("error");
    pl_free(pl);
    return EXIT_FAILURE;
  }
  printf("sk : ");
  for (int i = 0; i < n; i++)
    printf("%d ",sk[i]);
  puts("\n");
  pl_print(pl);

  distribution(l,t,pl);

  // choice of the r participants who want to recover the secret vector
  int tab[n], len = n, r = n - t;
  for (int i = 0; i < len; i++)
    tab[i] = i + 1;
  for (int i = 0; i < r; i++) {
    int ind = rand() % len;
    int in = tab[ind];
    pl->sigg[0][i] = in;
    pl->sigg[1][i] = mod((pl->sig[in])^inv(sk[in],q),q);
    pl->DLEQ[i] = 0;
    len--;
    for (int j = ind; j < len - 1; j++)
      tab[j] = tab[j+1];
  }

  printf("r = %d\n", r);
  for (int i = 0; i < r; i++)
    printf("%d ",pl->sigg[0][i]);


  reconstruction(r, pl);

  pl_print(pl);

  free(sk);
  pl_free(pl);


  // // test lambda
  // puts("test lambda");
  // int tab[5] = {1,2,3,4,5};
  // for (int i = 0; i < 5; i++) {
  //   for (int j = 0; j < 3; j++) {
  //     printf("lambda(%d,%d) = %d\n",tab[i],j,lambda(tab[i],j,tab,5));
  //   }
  // }

  return EXIT_SUCCESS;
}


/* TODO :
- régler problème de reconstruction
  -> voir si sigg correspondent à h^sig 
- implémentation figure 6 dans scrape++
- nettoyer le code */
