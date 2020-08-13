#include <iostream>
#include <cstdlib>
#include <cstdbool>
#include <vector>
#include <ctime>
#include <unistd.h>

#include "pvss.hpp"
#include "func.hpp"


bool dist = false;
bool rec = false;


////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// public ledger /////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

struct pl_t {
  int n; // number of participants
  int t; // threshold
  int l; // number of secrets
  int r; // number of participants that wanting to reconstruct the secrets
  ec_t* pk; // publi keys
  bn_t* sig; // shamir's shares
  ec_t* sighat; // encrypted shares
  int *reco_parties;
  ec_t* sigtilde; // decrypted shares and their index
  ec_t* S; // secrets reconstructed
};


pl_t *pl_alloc(const int n) {
  pl_t *pl;
  pl = new (nothrow) pl_t;
  if (!pl)
    return NULL;
  pl->n = n;
  pl->pk = new ec_t[n];
  return pl;
}

void pl_free(pl_t *pl) {
  if (pl) {
    delete[] pl->pk;
    if (dist) {
      delete[] pl->sig;
      delete[] pl->sighat;
    }
    if (rec) {
      delete[] pl->reco_parties;
      delete[] pl->sigtilde;
      delete[] pl->S;
    }
    delete pl;
  }
}

void pl_print(pl_t *pl) {
  cout << "\n\n___________________________ Public Ledger _________________________\n";
  cout << pl->n << " participants" << endl;
  cout << "pk :" << endl;
  for (int i = 0; i <pl->n; ++i)
    ec_print(pl->pk[i]);
  if (dist) {
    cout << endl << pl->t << " threshold" << endl;
    cout << "shamir's shares :" << endl;
    for (int i = 0; i < pl->n; ++i)
      bn_print(pl->sig[i]);
    cout << "\nencrypted shares :" << endl;
    for (int i = 0; i < pl->n; ++i)
      ec_print(pl->sighat[i]);
  }
  if (rec) {
    cout << endl << pl->r << " parties want to reconstruct" << endl;
    cout << "decrypted shares :" << endl;
    for (int i = 0; i < pl->r; i++) {
      cout << pl->reco_parties[i] << endl;
      ec_print(pl->sigtilde[i]);
    }
    cout << "\n\nThe " << pl->l << " secrets are :\n";
    for (int i = 0; i < pl->l; ++i) {
      cout << "S" << pl->l-i-1 << endl;
      ec_print(pl->S[i]);
    }
  }
  cout << "\n___________________________________________________________________\n\n";
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// scheme //////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

pl_t *setup(bn_t* sk, const int n, bn_t q) {
  for (int i = 0; i < n; i++) {
    bn_new(sk[i]);
    bn_rand_mod(sk[i],q);
    while (bn_is_zero(sk[i])) {
      bn_new(sk[i]);
      bn_rand_mod(sk[i],q);
    }
  }
  pl_t * pl = pl_alloc(n);
  if (!pl)
    return NULL;
  for (int i = 0; i < n; i++) {
    ec_mul_gen(pl->pk[i],sk[i]);
  }
  pl->n = n;
  bn_free(tmp);
  return pl;
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

void distribution(const int l, const int t, pl_t *pl, bn_t q) {
  if (!pl || t < 1 || t > pl->n)
    return;
  int deg = t + l;
  // random polynomila of Fp
  bn_t *P = new bn_t[deg];
  for (int i = 0; i < deg; ++i) {
    bn_rand_mod(P[i],q);
  }
  pl->l = l;
  pl->t = t;
  bn_t *s = new bn_t[pl->n+l];
  for (int i = -l+1; i <= pl->n; i++) {
    apply_poly(s[i+l-1],P,deg,i,q);
  }

  // attribution of the shamir's shares their values, computation of encrypted shares
  // and proof ldei
  pl->sig = new bn_t[pl->n];
  pl->sighat = new ec_t[pl->n];
  for (int i = 0; i < pl->n; i++) {
    bn_copy(pl->sig[i],s[i+l]);
    ec_mul(pl->sighat[i],pl->pk[i],s[i+l]);
  }
  // pl->ld.prove(pl->q, p, pl->pk, alpha, deg, pl->sighat, P);

  // print secrets for verification //////////////////////////
  // cout << endl << "Secrets :\n";
  // ec_t etmp;
  // ec_null(etmp);
  // for (int i = l-1; i >= 0; i--) {
  //   ec_new(etmp);
  //   ec_mul_gen(etmp, s[l-1-i]);
  //   cout << "S" << i << " =\n";
  //   ec_print(etmp);
  // }
  dist = true;
  // CLEAN UP
  delete[] s;
  delete[] P;
  ec_free(etmp);
}

void lambda(bn_t** lambs, const int t, pl_t *pl, bn_t q) {
  bn_t num, den, tmp1, tmp2;
  bn_null(num);
  bn_null(den);
  bn_null(tmp1);
  bn_null(tmp2);
  for (int j = 0; j < pl->l; j++) {
    bn_new(num);
    bn_new(den);
    for (int i = 0; i < t; i++) {
      bn_set_dig(num,1);
      bn_set_dig(den,1);
      for (int m = 0; m < t; m++)
        if (m != i) {
          bn_new(tmp1);
          bn_set_dig(tmp1, j);
          bn_neg(tmp1,tmp1);
          bn_new(tmp2);
          bn_sub_dig(tmp2,tmp1,pl->reco_parties[m]);
          bn_mul(num, num, tmp2);
          bn_mod(num,num,q);
          bn_new(tmp1);
          bn_set_dig(tmp1, pl->reco_parties[i]);
          bn_new(tmp2);
          bn_sub_dig(tmp2, tmp1, pl->reco_parties[m]);
          bn_mul(den, den, tmp2);
          bn_mod(den,den,q);
        }
      bn_new(tmp1);
      bn_mod_inv(tmp1,den,q);
      bn_mod(num,num,q);
      bn_new(tmp2);
      bn_mul(tmp2,num,tmp1);
      bn_mod(lambs[i][j],tmp2,q);
    }
  }
  bn_free(num);
  bn_free(den);
  bn_free(tmp1);
  bn_free(tmp2);
}

void reconstruction(const int r, pl_t *pl, bn_t q) {
  if (!pl) // error : the public ledger doesn't exist
    return;
  int t = pl->n - pl->t;
  if (r < t) // error : not enough parts
    return;
  pl->S = new ec_t[pl->l];
  bn_t** lambs;
  lambs = new bn_t*[t];
  for (int i = 0; i < t; ++i)
    lambs[i] = new bn_t[pl->l];
  clock_t time_lambda = clock();
  lambda(lambs,t,pl,q);
  time_lambda = clock() - time_lambda;
  ec_t tmp;
  ec_null(tmp);
  for (int j = 0; j < pl->l; j++) {
    // ec_new(tmp);
    // ec_mul(tmp,pl->sigtilde[0],lambs[0][j]);
    // ec_copy(pl->S[pl->l-j-1], tmp);
    ec_new(pl->S[pl->l-j-1]);
    ec_mul(pl->S[pl->l-j-1],pl->sigtilde[0],lambs[0][j]);
    for (int i = 1; i < t; i++) {
      ec_new(tmp);
      ec_mul(tmp,pl->sigtilde[i],lambs[i][j]);
      ec_add(pl->S[pl->l-j-1], pl->S[pl->l-j-1], tmp);
    }
  }
  for (int i = 0; i < t; ++i)
    delete[] lambs[i];
  delete[] lambs;
  pl->r = r;
  ec_free(tmp);
  rec = true;
  cout << "time for computing the lambdas: " << (float)time_lambda/CLOCKS_PER_SEC << "s" << endl;
}

void pvss_test(void) {

  /////////////////////////////////////////////////// initialization of Fq and E
  if (core_init() != RLC_OK) {
    core_clean();
    return;
  }

  // Ask for field to get Tweedledum.
	if (fp_param_set_any() == RLC_ERR) {
		RLC_THROW(ERR_NO_FIELD);
		core_clean();
		return;
	}

	// Ask for curve to get Tweedledum.
	if (ec_param_set_any() == RLC_ERR) {
		RLC_THROW(ERR_NO_CURVE);
		core_clean();
		return;
	}

  bn_t q;
  bn_new(q);
  ec_curve_get_ord(q);

  /////////////////////////////////////////////////////////// times declarations
  clock_t rec0, rec, setup_time, dist_time, decrypt_time, reco_time, all_time;
  rec0 = clock();

  /////////////////////////////////////////////////////////////////// PARAMETERS
  int n = 200;
  int t = n/3;
  int l = n-2*t;

  cout << "\n\ntimes for " << n << " participants with elliptic curves:\n\n";
  /////////////////////////////////////////////////////////////////////// SET UP
  bn_t* sk;
  sk = new bn_t[n];
  rec = clock();
  pl_t *pl = setup(sk,n,q);
  setup_time = clock() - rec;
  if (!pl)
    return;
  // pl_print(pl);


  ///////////////////////////////////////////////////////////////// DISTRIBUTION
  rec = clock();
  distribution(l,t,pl,q);
  dist_time = clock() - rec;

  ///////////////////////////////////// SHARE OF DECRYPTED SHARES AND PROOF DLEQ
  // choice of participant wanting to reconstruct
  int tab[n];
  int len = n;
  int r = n - t;
  bn_t tmp;
  bn_null(tmp);
  pl->sigtilde = new ec_t[r];
  pl->reco_parties = new int[r];
  bn_t* invsk = new bn_t[r];
  prng_init(time(NULL) + getpid());
  for (int i = 0; i < len; i++)
    tab[i] = i;
  for (int i = 0; i < r; i++) {
    int ind = rand() % len;
    int in = tab[ind];
    pl->reco_parties[i] = in + 1;
    bn_new(tmp);
    bn_mod_inv(tmp,sk[in],q);
    bn_copy(invsk[i],tmp);
    len--;
    for (int j = ind; j < len; j++)
      tab[j] = tab[j+1];
  }
  // computations
  rec = clock();
  int id;
  for (int i = 0; i < r; i++) {
    id = pl->reco_parties[i];
    ec_mul(pl->sigtilde[i],pl->sighat[id-1], invsk[i]);
  }
  decrypt_time = clock() - rec;

  /////////////////////////////////////////////////////////////// RECONSTRUCTION
  rec = clock();
  reconstruction(r,pl,q);
  reco_time = clock() - rec;
  // pl_print(pl);

  ///////////////////////////////////////////////////////////////////// CLEAN UP
  pl_free(pl);
  delete[] sk;
  delete[] invsk;
  core_clean();

  //////////////////////////////////////////////////////////////////////// TIMES
  all_time = clock() - rec0;
  cout << "time for seting up: " << (float)setup_time/CLOCKS_PER_SEC << "s" << endl;
  cout << "time for distributing: " << (float)dist_time/CLOCKS_PER_SEC << "s" << endl;
  cout << "time for sharing decrypted shares: " << (float)decrypt_time/CLOCKS_PER_SEC << "s" << endl;
  cout << "time for reconstructing the secrets: " << (float)reco_time/CLOCKS_PER_SEC << "s" << endl;
  cout << "\nglobal time: " << (float)all_time/CLOCKS_PER_SEC << "s" << endl << endl;
}
