/* implementation of pvss protocol */
#include <iostream>
#include <cstdlib>
#include <cstdbool>
#include <ctime>
#include <cmath>
#include <unistd.h>

#include <gmp.h>
#include <NTL/ZZ_pX.h>

#include "pvss.hpp"
#include "ldei.hpp"

// structures and variables
bool dist = false;
bool rec = false;

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
  ldei_t *LDEI; //proof LDEI
  Mat<ZZ_p> sigtilde; // decrypted shares and their index
  // int *DLEQ;
  Vec<ZZ_p> S; // secrets reconstructed
};

void prng_init(const unsigned int seed) {
  static bool seed_init = false;
  if (! seed_init) {
    srand(seed);
    seed_init = true;
  }
}

pl_t *pl_alloc(const int n) {
  pl_t *pl;
  pl = new (nothrow) pl_t;
  if (!pl)
    return NULL;
  pl->n = n;
  pl->t = 0;
  pl->q = ZZ(0);
  pl->r = 0;
  pl->l = 0;
  pl->pk.SetLength(n);
  pl->sig.SetLength(n);
  pl->sighat.SetLength(n);
  pl->LDEI = ldei_alloc(n);
  return pl;
}

void pl_free(pl_t *pl) {
  if (pl) {
    pl->q.kill();
    pl->pk.kill();
    pl->sig.kill();
    pl->sighat.kill();
    pl->sigtilde.kill();
    ldei_free(pl->LDEI);
    delete pl;
  }
}

void pl_print(pl_t *pl) {
  cout << "\n\n___________________________ Public Ledger ___________________________\n";
  cout << pl->n << " participants" << endl;
  cout << "q = " << pl->q << endl;
  cout << "pk : " << pl->pk << endl;
  if (dist) {
    cout << endl << pl->t << " threshold" << endl;
    cout << "shamir's shares :" << endl << pl->sig << endl;
    cout << "encrypted shares :" << endl << pl->sighat << endl;
    if (pl->LDEI)
      ldei_print(pl->LDEI);
    cout << endl;
  }
  if (rec) {
    cout << endl << pl->r << " parties want to reconstruct" << endl;
    cout << "decrypted shares :" << endl;
    for (int i = 1; i <= pl->r; i++)
      cout << pl->sigtilde(i) << " ";
    // if (pl->DLEQ) {
    //   cout << endl << "      DLEQ : ";
    //   for (int i = 0; i < pl->r; i++)
    //     cout << pl->DLEQ[i] << " ";
    // }
    cout << "\n\nThe " << pl->l << " secrets are :\n" << pl->S << endl;
  }
  cout << "_____________________________________________________________________" << endl;
}

// in ZZp
void generator(ZZ_p& g, const ZZ& q) {
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

pl_t *setup(Vec<ZZ_p>& sk, const int n, const ZZ& q, const ZZ& p, const ZZ_p& h) {
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
  pl_t * pl = pl_alloc(n);
  if (!pl)
      return NULL;
  for (int i = 0; i < n; i++) {
    r = rep(sk[i]);
    power(pl->pk[i], h, r);
  }
  pl->n = n;
  pl->h = h;
  pl->q = q;
  return pl;
}

void distribution(const int l, const int t, const Vec<ZZ_p>& alpha, pl_t *pl) {
  if (!pl || t < 1 || t > pl->n)
    return;
  int deg = t + l;
  ZZ_pX P;
  random(P, deg);
  pl->l = l;
  pl->t = t;
  Vec<ZZ_p> s;
  s.SetLength(pl->n+l);
  ZZ_p tmp;
  ZZ repzz;
  for (int i = -l+1; i <= pl->n; i++) {
    tmp = ZZ_p(i);
    eval(s[i+l-1], P, tmp);
  }
  ZZ p = 2 * pl->q + 1;
  ZZ_pPush push(p);
  // attribution of the shamir's shares their values, computation of encrypted shares
  // and proof ldei
  for (int i = 0; i < pl->n; i++) {
    pl->sig[i] = s[i+l];
    repzz = rep(s[i+l]);
    power(pl->sighat[i],pl->pk[i],repzz);
  }
  pl->LDEI = ldei_prove(pl->q, p, pl->pk, alpha, deg, pl->sighat, P);
  // print secrets for verification //////////////////////////
  // cout << endl << "Secrets :";
  // computation of secrets
  for (int i = l-1; i >= 0; i--) {
    repzz = rep(s[i]);
    power(tmp,pl->h,repzz);
    // cout << "S" << i << " = h^" << repzz << " = " << tmp << endl;
  }
  dist = true;
  s.kill();
}

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

void reconstruction(const int r, pl_t *pl) {
  if (!pl) // error : the public ledger doesn't exist
    return;
  int t = pl->n - pl->t;
  if (r < t) // error : not enough parts
    return;
  // verification of proof DLEQ
  pl->S.FixLength(pl->l);
  ZZ_pPush push(pl->q);
  Mat<ZZ_p> lambs;
  lambda(lambs, t, pl);
  ZZ lamb;
  ZZ_p::init(2 * pl->q + 1);
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
  lambs.kill();
}

void pvss(void) {

  clock_t rec0, rec, setup_time, dist_time, ldeiverif_time, decrypt_time,
    dleqverif_time, reco_time, recoverif_time, all_time;
  rec0 = clock();

  // PARAMETERS
  int n = 1024;
  ZZ q;
  GenGermainPrime(q,1024);
  ZZ p = 2 * q + 1;
  int t = n/3;
  int l = n-2*t;
  ZZ_p::init(p);
  ZZ_p g;
  generator(g,p);
  ZZ_p h;
  power(h,g,2);

  // SET UP
  ZZ_p::init(q);
  Vec<ZZ_p> sk;
  rec = clock();
  pl_t *pl = setup(sk,n,q,p,h);
  setup_time = clock() - rec;
  if (!pl)
    return;
  // cout << endl << "sk : " << sk << endl;


  // DISTRIBUTION
  Vec<ZZ_p> alpha;
  alpha.SetLength(pl->n);
  for (int i = 0; i < pl->n; i++)
    alpha[i] = ZZ_p(i+1);
  rec = clock();
  distribution(l,t,alpha,pl);
  dist_time = clock() - rec;

  // VERIFICATION
  rec = clock();
  if (ldei_verify(pl->LDEI, q, p, pl->pk, alpha, t+l, pl->sighat) == false)
    return;
  ldeiverif_time = clock() - rec;

  // SHARE OF DECRYPTED SHARES AND PROOF DLEQ
  // choice of participant wanting to reconstruct
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
    // cout << in << " ";
    inv(invsk[i-1],sk[in]);
    choice[i] = in;
    len--;
    for (int j = ind; j < len; j++)
      tab[j] = tab[j+1];
  }
  // computations
  ZZ_p::init(p);
  rec = clock();
  for (int i = 1; i <= r; i++) {
    power(tmp,pl->sighat[choice[i]],rep(invsk[i-1]));
    pl->sigtilde(i,2) = tmp;
    // pl->DLEQ[i] = 0;
  }
  decrypt_time = clock() - rec;
  invsk.kill();


  // VERIFICATION OF PROOF DLEQ
  rec = clock();
  dleqverif_time = clock() - rec;

  // RECONSTRUCTION
  rec = clock();
  reconstruction(r, pl);
  reco_time = clock() - rec;
  pl_print(pl);

  // RECONSTRUCTINO VERIFICATION
  rec = clock();
  recoverif_time = clock() - rec;


  // CLEAN UP
  pl_free(pl);
  q.kill();

  // TIME
  all_time = clock() - rec0;
  cout << "time for seting up: " << (float)setup_time/CLOCKS_PER_SEC << "s" << endl;
  cout << "time for distributing: " << (float)dist_time/CLOCKS_PER_SEC << "s" << endl;
  cout << "time for verifying ldei: " << (float)ldeiverif_time/CLOCKS_PER_SEC << "s" << endl;
  cout << "time for sharing decrypted shares with dleq: " << (float)decrypt_time/CLOCKS_PER_SEC << "s" << endl;
  cout << "time for verifying dleq: " << (float)dleqverif_time/CLOCKS_PER_SEC << "s" << endl;
  cout << "time for reconstructing the secrets: " << (float)reco_time/CLOCKS_PER_SEC << "s" << endl;
  cout << "time for verifying the reconstruction: " << (float)recoverif_time/CLOCKS_PER_SEC << "s" << endl;
  cout << "time for doing all these things: " << (float)all_time/CLOCKS_PER_SEC << "s" << endl;
}
