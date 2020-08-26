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

clock_t setup(pl_t* pl, bn_t* sk, const int n, bn_t q) {
  clock_t time = 0, timetmp;
  // choice of secret keys sk
  for (int i = 0; i < n; i++) {
    bn_new(sk[i]);
    bn_rand_mod(sk[i],q);
    while (bn_is_zero(sk[i])) {
      bn_new(sk[i]);
      bn_rand_mod(sk[i],q);
    }
  }
  // computation of public keys pk = h^sk
  for (int i = 0; i < n; i++) {
    timetmp = clock();
    ec_mul_gen(pl->pk[i],sk[i]);
    time += clock() - timetmp;
  }
  // set the variables in the public ledger
  pl->n = n;
  bn_free(tmp);
  return time;
}

clock_t distribution(const int l, const int t, pl_t *pl, bn_t q) {
  if (!pl || t < 1 || t > pl->n)
    return 0;
  int deg = t + l;
  // choice of the polynomial P
  bn_t *P = new bn_t[deg];
  for (int i = 0; i < deg; ++i) {
    bn_rand_mod(P[i],q);
  }
  // conputation of the exponents and shamir's shares
  bn_t *s = new bn_t[pl->n+l];
  for (int i = -l+1; i <= pl->n; i++) {
    apply_poly(s[i+l-1],P,deg,i,q);
  }

  clock_t time = 0, timetmp;
  // computation of encrypted shares
  pl->sighat = new ec_t[pl->n];
  for (int i = 0; i < pl->n; i++) {
    timetmp = clock();
    ec_mul(pl->sighat[i],pl->pk[i],s[i+l]);
    time += clock() - timetmp;
  }
  // set the variables in the public ledger
  pl->l = l;
  pl->t = t;
  dist = true;
  // clean up
  delete[] s;
  delete[] P;
  ec_free(etmp);
  return time;
}

clock_t reconstruction(const int r, pl_t *pl, bn_t q) {
  if (!pl) // error : the public ledger doesn't exist
    return 0;
  int t = pl->n - pl->t;
  if (r < t) // error : not enough parts
    return 0;
  // computation of the Lagrange coefficients
  bn_t** lambs;
  lambs = new bn_t*[t];
  for (int i = 0; i < t; ++i)
    lambs[i] = new bn_t[pl->l];
  clock_t time = 0, timetmp;
  lambda(lambs,t,pl,q);
  // copmputation of the secrets
  pl->S = new ec_t[pl->l];
  ec_t tmp;
  ec_null(tmp);
  for (int j = 0; j < pl->l; j++) {
    ec_new(pl->S[pl->l-j-1]);
    timetmp = clock();
    ec_mul(pl->S[pl->l-j-1],pl->sigtilde[0],lambs[0][j]);
    time += clock() - timetmp;
    for (int i = 1; i < t; i++) {
      ec_new(tmp);
      timetmp = clock();
      ec_mul(tmp,pl->sigtilde[i],lambs[i][j]);
      ec_add(pl->S[pl->l-j-1], pl->S[pl->l-j-1], tmp);
      time += clock() - timetmp;
    }
  }
  // set the variables in the public ledger
  pl->r = r;
  rec = true;
  // clean up
  for (int i = 0; i < t; ++i)
    delete[] lambs[i];
  delete[] lambs;
  ec_free(tmp);
  return time;
}


////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// ppvss test ////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void pvss_test(const int n) {

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
  clock_t timetmp, setup_time, dist_time, decrypt_time = 0, reco_time, all_time;

  /////////////////////////////////////////////////////////////////// PARAMETERS
  int t = n/3;
  int l = n-2*t;

  /////////////////////////////////////////////////////////////////////// SET UP
  bn_t* sk;
  sk = new bn_t[n];
  pl_t * pl = pl_alloc(n);
  if (!pl)
    return;
  setup_time = setup(pl,sk,n,q);
  if (!pl)
    return;
  // pl_print(pl);


  ///////////////////////////////////////////////////////////////// DISTRIBUTION
  dist_time = distribution(l,t,pl,q);

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
  int id;
  for (int i = 0; i < r; i++) {
    id = pl->reco_parties[i];
    timetmp = clock();
    ec_mul(pl->sigtilde[i],pl->sighat[id-1], invsk[i]);
    decrypt_time += clock() - timetmp;
  }

  /////////////////////////////////////////////////////////////// RECONSTRUCTION
  reco_time = reconstruction(r,pl,q);
  // pl_print(pl);

  ///////////////////////////////////////////////////////////////////// CLEAN UP
  pl_free(pl);
  delete[] sk;
  delete[] invsk;
  core_clean();

  all_time = setup_time + dist_time + decrypt_time + reco_time;
  //////////////////////////////////////////////////////////////////////// TIMES
  cout << "\n\ntimes for " << n << " participants with elliptic curves:\n\n";
  cout << "time for setting up: " << (float)setup_time/CLOCKS_PER_SEC << "s\n";
  cout << "time for distributing: " << (float)dist_time/CLOCKS_PER_SEC << "s\n";
  cout << "time for sharing decrypted shares: " << (float)decrypt_time/CLOCKS_PER_SEC << "s\n";
  cout << "time for reconstructing the secrets: " << (float)reco_time/CLOCKS_PER_SEC << "s\n";
  cout << "\nglobal time: " << (float)all_time/CLOCKS_PER_SEC << "s\n\n";
  cout << "here we only timed elementary operations and we didn't time the proofs\n\n";
}
