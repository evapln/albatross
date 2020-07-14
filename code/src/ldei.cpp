#include <cstdbool>

#include "ldei.hpp"
#include "hash.hpp"

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// LDEI ////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

struct ldei_t {
  Vec<ZZ_p> a;
  ZZ_p e;
  ZZ_pX z;
};

ldei_t *ldei_alloc(const int m) {
  ldei_t *ld;
  ld = new (nothrow) ldei_t;
  if (!ld)
    return NULL;
  ld->a.SetLength(m);
  ld->e = ZZ_p(0);
  ld->z = 0;
  return ld;
}

/* v only for test v */
ldei_t* ldei_copy(ldei_t* src, int m) {
  ldei_t* dest = ldei_alloc(m);
  dest->a = src->a;
  dest->e = src->e;
  dest->z = src->z;
  return dest;
}

void ldei_set_a(ldei_t* ld, Vec<ZZ_p>& a) {
  ld->a = a;
}

void ldei_set_e(ldei_t* ld, ZZ_p& e) {
  ld->e = e;
}

void ldei_set_z(ldei_t* ld, ZZ_pX& z) {
  ld->z = z;
}
/* ^ only for test ^ */

void ldei_free(ldei_t *ld) {
  if (ld) {
    ld->a.kill();
    delete ld;
  }
}



void ldei_print(ldei_t *ld) {
  cout << "\n a = " << ld->a << "\n e = " << ld->e << "\n z = " << ld->z << endl;
}

// in ZZq
ldei_t* ldei_prove(const ZZ& q, const ZZ& p, const Vec<ZZ_p>& g, const Vec<ZZ_p>& alpha,
const long k, const Vec<ZZ_p>& x, const ZZ_pX& P) {
  ZZ_p::init(q);
  // verify entries
  int m = g.length();
  if (alpha.length() != m || x.length() != m || k > m)
    return NULL;


  //create the vector a
  ZZ_pX R = random_ZZ_pX(k+1);
  Vec<ZZ_p> R_eval = eval(R, alpha);

  ZZ_p::init(p);

  // create ldei structure
  ldei_t *ld = ldei_alloc(m);

  for (int i = 0; i < m; i++) {
    power(ld->a[i], g[i], rep(R_eval[i]));
  }

  // in ZZq
  ZZ_pPush push(q);

  // create e the digest of the hash function
  hash_ZZp(ld->e, x, ld->a);

  // create z the polynomial of ld
  ZZ_pX tmp = ld->e * P;
  ld->z = tmp + R;

  return ld;
}

bool ldei_verify(const ldei_t* ld, const ZZ& q, const ZZ& p, const Vec<ZZ_p>& g,
const Vec<ZZ_p>& alpha, const long k, const Vec<ZZ_p>& x){
  // test of entries
  if (!ld)
    return false;
  int m = ld->a.length();
  if (alpha.length() != m || x.length() != m || g.length() != m) {
    cout << "Verification not passed: bad length\n";
    return false;
  }

  // test degree polynomial
  if (deg(ld->z) > k) {
    cout << "Verification not passed: bad degree of z\n";
    return false;
  }

  ZZ_p::init(q);

  // test of e = hash(x,a)
  ZZ_p hashzz;
  hash_ZZp(hashzz, x, ld->a);
  if (ld->e != hashzz) {
    cout << "Verification not passed: bad digest\n";
    return false;
  }


  // test of x_i^e * a_i = g_i^z(alpha_i)
  Vec<ZZ_p> zi = eval(ld->z, alpha);
  ZZ_pPush push(p);
  ZZ_p tmp4, tmp2, tmp3;
  for (int i = 0; i < m; i++) {
    power(tmp2, g[i], rep(zi[i]));
    power(tmp3, x[i], rep(ld->e));
    mul(tmp4, tmp3, ld->a[i]);
    if (tmp2 != tmp4) {
      cout << "Verification not passed: bad a_i\n";
      return false;
    }
  }
  cout << "Verification passed\n";
  return true;
}


////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// local LDEI ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// in ZZq
bool localldei(const Vec<ZZ_p>& a, const long k, const Vec<ZZ_p>& x,
const long m, const ZZ& p) {

  // computation of vector u
  Vec<ZZ_p> u;
  u.SetLength(m);
  ZZ_p prod, tmp;
  for (int i = 0; i < m; i++) {
    prod = 1;
    for (int l = 0; l < m; l++)
      if (l != i) {
        sub(tmp, a[i], a[l]);
        prod *= tmp;
      }
    inv(u[i],prod);
  }

  // random polynomial
  ZZ_pX P;
  do random(P,m-k-1); while(IsZero(P));

  // computation of v
  Vec<ZZ_p> v;
  v.SetLength(m);
  for (int i = 0; i < m; i++) {
    eval(tmp, P, a[i]);
    mul(v[i], u[i], tmp);
  }



  // verification
  ZZ_pPush push(p);
  prod = 1;
  for (int i = 0; i < m; i++) {
    power(tmp, x[i], rep(v[i]));
    prod *= tmp;
  }

  return (IsOne(prod));
}
