#include <iostream>
#include <cstdlib>
#include <cstdbool>
#include <vector>
#include <gmpxx.h>
extern "C" {
  #include "relic.h"
}

bool dist = false;
bool rec = false;

using namespace std;

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// public ledger /////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

struct pl_t {
  int n; // number of participants
  int t; // threshold
  int l; // number of secrets
  int r; // number of participants that wanting to reconstruct the secrets
  // ZZ q; // order of the group G_q
  // ec_t h; // generator of G_q
  ec_t* pk; // publi keys
  fp_t* sig; // shamir's shares
  ec_t* sighat; // encrypted shares
  // LDEI ld; //proof LDEI
  int *reco_parties;
  ec_t* sigtilde; // decrypted shares and their index
  // vector<DLEQ> dl; // proof DLEQ
  ec_t* S; // secrets reconstructed
};


pl_t *pl_alloc(const int n) {
  pl_t *pl;
  pl = new (nothrow) pl_t;
  if (!pl)
    return NULL;
  pl->n = n;
  pl->pk = new ec_t[n];
  // pl->sig = new fp_t[n];
  // pl->sighat = new ec_t[n];
  // pl->sigtilde.reserve(n);
  return pl;
}

void pl_free(pl_t *pl) {
  if (pl) {
    // for (int i = 0; i < pl->n; ++i) {
      // ec_free(pl->pk[i]);
    // }
    delete pl->pk;
    if (dist) {
      // for (int i = 0; i < pl->n; ++i) {
      //   fp_free(pl->sig[i]);
      //   ec_free(pl->sighat[i]);
      // }
      delete pl->sig;
      delete pl->sighat;
    }
    if (rec) {
      // for (int i = 0; i < pl->r; ++i) {
      //   ec_free(pl->sigtilde[i]);
      // }
      delete pl->reco_parties;
      delete pl->sigtilde;
      // for (int i = 0; i < pl->l; ++i) {
      //   ec_free(pl->S[i]);
      // }
      delete pl->S;
    }
    delete pl;
  }
}

void pl_print(pl_t *pl) {
  cout << "\n\n___________________________ Public Ledger _________________________\n";
  cout << pl->n << " participants" << endl;
  // cout << "q = " << pl->q << endl;
  cout << "pk :" << endl;
  for (int i = 0; i <pl->n; ++i)
    ec_print(pl->pk[i]);
  if (dist) {
    cout << endl << pl->t << " threshold" << endl;
    cout << "shamir's shares :" << endl;
    for (int i = 0; i < pl->n; ++i)
      fp_print(pl->sig[i]);
    cout << "\nencrypted shares :" << endl;
    for (int i = 0; i < pl->n; ++i)
      ec_print(pl->sighat[i]);
    // pl->ld.print();
    // cout << endl;
  }
  if (rec) {
    cout << endl << pl->r << " parties want to reconstruct" << endl;
    cout << "decrypted shares :" << endl;
    for (int i = 0; i < pl->r; i++) {
      cout << pl->reco_parties[i];
      ec_print(pl->sigtilde[i]);
    }
    // cout << endl << pl->sigtilde << endl;
    // int len = pl->dl.size();
    // for (int i = 0; i < len; i++) {
    //   cout << "dleq number " << i << endl;
    //   pl->dl[i].print();
    // }
    cout << "\n\nThe " << pl->l << " secrets are :\n";
    for (int i = 0; i < pl->l; ++i)
      ec_print(pl->S[i]);
  }
  cout << "\n___________________________________________________________________\n\n";
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// scheme //////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

pl_t *setup(fp_t* sk, const int n) {
  for (int i = 0; i < n; i++) {
    fp_rand(sk[i]);
    while (fp_is_zero(sk[i]))
      fp_rand(sk[i]);
  }
  pl_t * pl = pl_alloc(n);
  bn_t tmp;
  bn_null(tmp);
  if (!pl)
    return NULL;
  for (int i = 0; i < n; i++) {
    bn_new(tmp);
    fp_prime_back(tmp,sk[i]);
    ec_mul_gen(pl->pk[i],tmp);
  }
  pl->n = n;
  bn_free(tmp);
  return pl;
}

// return
void apply_poly(fp_t y, const fp_t *P, const int deg, const int x) {
  fp_t pow, tmp;
  fp_new(y);
  fp_new(pow);
  fp_new(tmp);
  fp_copy(y, P[0]);
  fp_set_dig(pow, x);
  for (int i = 1; i < deg; ++i) {
    fp_mul(tmp, P[i], pow);
    fp_add(y,y,tmp);
    fp_mul_dig(pow,pow,x);
  }
  fp_free(pow);
  fp_free(tmp);
}

void distribution(const int l, const int t, pl_t *pl) {
  if (!pl || t < 1 || t > pl->n)
    return;
  int deg = t + l;
  // random polynomila of Fp
  fp_t *P = new fp_t[deg];
  for (int i = 0; i < deg; ++i)
    fp_rand(P[i]);
  pl->l = l;
  pl->t = t;
  fp_t *s = new fp_t[pl->n+l];
  for (int i = -l+1; i <= pl->n; i++) {
    apply_poly(s[i+l-1], P, deg, i);
  }
  // attribution of the shamir's shares their values, computation of encrypted shares
  // and proof ldei
  bn_t tmp;
  bn_null(tmp);
  pl->sig = new fp_t[pl->n];
  pl->sighat = new ec_t[pl->n];
  for (int i = 0; i < pl->n; i++) {
    fp_copy(pl->sig[i],s[i+l]);
    bn_new(tmp);
    fp_prime_back(tmp,s[i+l]);
    ec_mul(pl->sighat[i],pl->pk[i],tmp);
  }
  // pl->ld.prove(pl->q, p, pl->pk, alpha, deg, pl->sighat, P);

  // print secrets for verification //////////////////////////
  cout << endl << "Secrets :\n";
  ec_t etmp;
  ec_null(etmp);
  // computation of secrets
  for (int i = l-1; i >= 0; i--) {
    ec_new(etmp);
    bn_new(tmp);
    fp_prime_back(tmp,s[l-1-i]);
    ec_mul_gen(etmp, tmp);
    cout << "S" << i << " =\n";
    ec_print(etmp);
  }
  dist = true;
  // CLEAN UP
  delete s;
  delete P;
  bn_free(tmp);
  ec_free(etmp);
}

void lambda(fp_t** lambs, const int t, pl_t *pl) {
  lambs = new fp_t*[t];
  for (int i = 0; i < t; ++i)
    lambs[i] = new fp_t[pl->l];
  fp_t num, den, tmp1, tmp2;
  fp_new(num);
  fp_new(den);
  fp_new(tmp1);
  fp_new(tmp2);
  for (int j = 0; j < pl->l; j++) {
    for (int i = 0; i < t; i++) {
      fp_prime_conv_dig(num,1);
      fp_prime_conv_dig(den,1);
      for (int m = 0; m < t; m++)
        if (m != i) {
          fp_set_dig(tmp1, -j);
          fp_sub_dig(tmp2,tmp1,pl->reco_parties[m]);
          fp_mul(num, num, tmp2);
          fp_set_dig(tmp1, pl->reco_parties[i]);
          fp_sub_dig(tmp2, tmp1, pl->reco_parties[m]);
          fp_mul(den, den, tmp2);
        }
      fp_inv(tmp1,den);
      fp_mul(lambs[i][j],num,tmp1);
    }
  }
  fp_free(num);
  fp_free(den);
  fp_free(tmp1);
  fp_free(tmp2);
}

void reconstruction(const int r, pl_t *pl) {
  if (!pl) // error : the public ledger doesn't exist
    return;
  int t = pl->n - pl->t;
  if (r < t) // error : not enough parts
    return;
  pl->S = new ec_t[pl->l];
  fp_t** lambs;
  lambda(lambs, t, pl);
  // fp_t lamb;
  ec_t tmp;
  bn_t lamb;
  bn_new(lamb);
  // fp_new(lamd);
  ec_new(tmp);
  for (int j = 0; j < pl->l; j++) {
    fp_prime_back(lamb, lambs[0][j]);
    ec_mul(tmp,pl->sigtilde[0],lamb);
    ec_copy(pl->S[pl->l-j-1], tmp);
    for (int i = 1; i < t; i++) {
      fp_prime_back(lamb, lambs[i][j]);
      ec_mul(tmp,pl->sigtilde[i],lamb);
      ec_add(pl->S[pl->l-j-1], pl->S[pl->l-j-1], tmp);
    }
  }
  pl->r = r;
  rec = true;
  delete lambs;
}


int main(void) {

	if (core_init() != RLC_OK) {
		core_clean();
		return 1;
	}

	util_banner("Tests for the EC module:", 0);

	/* Ask for field to get Tweedledum. */
	if (fp_param_set_any() == RLC_ERR) {
		RLC_THROW(ERR_NO_FIELD);
		core_clean();
		return 0;
	}
	fp_param_print();

	/* Ask for curve to get Tweedledum. */
	if (ec_param_set_any() == RLC_ERR) {
		RLC_THROW(ERR_NO_CURVE);
		core_clean();
		return 0;
	}

	/* Print curve name to make sure it is the right one. */
	ec_param_print();

  int n = 10, t = 3;
  fp_t* sk;
  sk = new fp_t[n];
  pl_t* pl = setup(sk,n);
  pl_print(pl);
  for (int i = 0; i < n; ++i)
    fp_print(sk[i]);

  distribution(n,t,pl);

  pl_print(pl);

  pl_free(pl);
  delete sk;

	// /* Check that generator has the right order. */
	// ec_curve_get_gen(g);
	// ec_curve_get_ord(r);
	// ec_mul(g, g, r);
	// // assert(ec_is_infty(g));
  //
	// /* Generate some random point and exponents to illustrate scalar mult. */
	// ec_rand(a);
	// ec_rand(b);
	// bn_rand_mod(l, r);
	// ec_mul(a, a, l);
	// ec_print(a);
  //
	// // util_banner("All tests have passed.\n", 0);
  //
  //
  // ec_t* points;
  // points = new ec_t[10];
  // for (int i = 0; i <10; ++i) {
  //   ec_new(points[i]);
  //   ec_rand(points[i]);
  //   // points.push_back(pi);
  // }
  // cout << "points" << endl;
  // for (int i = 0; i < 10; ++i) {
  //   ec_print(points[i]);
  //   // ec_free(points[i]);
  // }
  //
  // ec_t gen;
  // ec_new(gen);
  // ec_curve_get_gen(g);
  // cout << "gen" << endl;
  // ec_print(g);
  //
  // delete points;
  // ec_free(pi);
	// ec_free(g);
	// ec_free(a);
	// ec_free(b);
	// bn_free(l);
	// bn_free(r);
	core_clean();
	return 0;
}
