#include <cstdbool>

#include "proofs.hpp"
#include "hash.hpp"

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// LDEI  ///////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

LDEI::LDEI() {
  a.SetLength(0);
}

void LDEI::print() {
  if (a.length() > 0)
    cout << "\n a = " << a << "\n e = " << e << "\n z = " << z << endl;
  else
    cout << "PAS DE PREUVE\n";
}

void LDEI::prove(const ZZ& q, const ZZ& p, const Vec<ZZ_p>& g, const Vec<ZZ_p>& alpha,
const long k, const Vec<ZZ_p>& x, const ZZ_pX& P) {
  ZZ_p::init(q);
  // verify entries
  int m = g.length();
  if (alpha.length() != m || x.length() != m || k > m) {
    a.SetLength(0);
  }
  else {
    //create the vector a
    ZZ_pX R = random_ZZ_pX(k+1);
    Vec<ZZ_p> R_eval = eval(R, alpha);
    ZZ_p::init(p);
    // create proof structure
    a.SetLength(m);
    for (int i = 0; i < m; i++)
      power(a[i], g[i], rep(R_eval[i]));
    // in ZZq
    ZZ_pPush push(q);
    // create e the digest of the hash function
    hash_ZZp(e, x, a);
    // create z the polynomial of ld
    ZZ_pX tmp = e * P;
    z = tmp + R;
    //clean up
    R_eval.kill();
  }
}

bool LDEI::verify(const ZZ& q, const ZZ& p, const Vec<ZZ_p>& g, const Vec<ZZ_p>& alpha,
const long k, const Vec<ZZ_p>& x) {
  // test of entries
  int m = a.length();
  if (alpha.length() != m || x.length() != m || g.length() != m) {
    cout << "Verification not passed: bad length\n";
    return false;
  }
  // test degree polynomial
  if (deg(z) > k) {
    cout << "Verification not passed: bad degree of z\n";
    return false;
  }
  ZZ_p::init(q);
  // test of e = hash(x,a)
  ZZ_p hashzz;
  hash_ZZp(hashzz, x, a);
  if (e != hashzz) {
    cout << "Verification not passed: bad digest\n";
    return false;
  }
  // test of x_i^e * a_i = g_i^z(alpha_i)
  Vec<ZZ_p> zi = eval(z, alpha);
  ZZ_pPush push(p);
  ZZ_p tmp1, tmp2, tmp3;
  for (int i = 0; i < m; i++) {
    power(tmp2, g[i], rep(zi[i]));
    power(tmp3, x[i], rep(e));
    mul(tmp1, tmp3, a[i]);
    if (tmp2 != tmp1) {
      cout << "Verification not passed: bad a_i\n";
      zi.kill();
      return false;
    }
  }
  zi.kill();
  return true;
}

////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// DLEQ ////// //////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

DLEQ::DLEQ() {
  a.SetLength(0);
}

void DLEQ::print() {
  if (a.length() > 0)
    cout << "\n a = " << a << "\n e = " << e << "\n z = " << z << endl;
  else
    cout << "PAS DE PREUVE\n";
}

void DLEQ::prove(const ZZ& q, const ZZ& p, const Vec<ZZ_p>& g, const Vec<ZZ_p>& x,
ZZ_p& alpha) {
  int m = g.length();
  if (x.length() != m) {
    a.SetLength(0);
  }
  else {
    ZZ_p::init(q);
    ZZ_p w = random_ZZ_p();
    ZZ_p::init(p);
    a.SetLength(m);
    for (int i = 0; i < m; i++) {
      power(a[i],g[i],rep(w));
    }
    ZZ_pPush push(q);
    hash_ZZp(e,g,x,a);
    ZZ tmp;
    mul(tmp,rep(alpha),rep(e));
    sub(z,rep(w),tmp);
  }
}

bool DLEQ::verify(const ZZ& q, const ZZ& p, const Vec<ZZ_p>& g, const Vec<ZZ_p>& x) {
  // test of entries
  int m = a.length();
  if (x.length() != m || g.length() != m) {
    cout << "Verification not passed: bad length\n";
    return false;
  }

  ZZ_p::init(q);

  // test of digest e
  ZZ_p hashzz;
  hash_ZZp(hashzz,g,x,a);
  if (e != hashzz) {
    cout << "Verification not passed: bad digest\n";
    return false;
  }

  // test of values of a_i = g_i^z * x_i^e
  ZZ_pPush push(p);
  ZZ_p tmp1, tmp2, tmp3;
  for (int i = 0; i < m; i++) {
    power(tmp2, g[i], z);
    power(tmp3, x[i], rep(e));
    mul(tmp1, tmp2, tmp3);
    if (a[i] != tmp1) {
      cout << "Verification not passed: bad a_i\n";
      return false;
    }
  }

  return true;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////// proof strucure //////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// struct proof_t {
//   Vec<ZZ_p> a;
//   ZZ_p e;
//   ZZ_pX z;
// };

// proof_t *proof_alloc(const int m) {
//   proof_t *ld;
//   ld = new (nothrow) proof_t;
//   if (!ld)
//     return NULL;
//   ld->a.SetLength(m);
//   ld->e = ZZ_p(0);
//   ld->z = 0;
//   return ld;
// }
//
// /* v only for test v */
// proof_t* proof_copy(proof_t* src, int m) {
//   proof_t* dest = proof_alloc(m);
//   dest->a = src->a;
//   dest->e = src->e;
//   dest->z = src->z;
//   return dest;
// }
//
// void proof_set_a(proof_t* ld, Vec<ZZ_p>& a) {
//   ld->a = a;
// }
//
// void proof_set_e(proof_t* ld, ZZ_p& e) {
//   ld->e = e;
// }
//
// void proof_set_z(proof_t* ld, ZZ_pX& z) {
//   ld->z = z;
// }
// /* ^ only for test ^ */
//
// void proof_free(proof_t *ld) {
//   if (ld) {
//     ld->a.kill();
//     delete ld;
//   }
// }

// void proof_print(proof_t *ld) {
  // cout << "\n a = " << ld->a << "\n e = " << ld->e << "\n z = " << ld->z << endl;
// }

// void Proof::print() {
//   cout << "\n a = " << ld->a << "\n e = " << ld->e << "\n z = " << ld->z << endl;
// }

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////// DLEQ //////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// proof_t* dleq_prove(const ZZ& q, const ZZ& p, const Vec<ZZ_p>& g,
// const Vec<ZZ_p>& x, ZZ_p& alpha) {
//   int m = g.length();
//   if (x.length() != m)
//     return NULL;
//   ZZ_p::init(q);
//   proof_t* dleq = proof_alloc(m);
//   ZZ_p w = random_ZZ_p();
//   ZZ_p::init(p);
//   for (int i = 0; i < m; i++)
//     power(dleq->a[i],g[i],rep(w));
//   ZZ_pPush push(q);
//   hash_ZZp(dleq->e,g,x,dleq->a);
//   ZZ_p tmp;
//   mul(tmp,alpha,dleq->e);
//   sub(dleq->z,w,tmp);
//   return dleq;
// }
//
// bool dleq_verify(const proof_t* dleq, const ZZ& q, const ZZ& p,
// const Vec<ZZ_p>& g, const Vec<ZZ_p>& x) {
  // // test of entries
  // if (!dleq)
  //   return false;
  // int m = dleq->a.length();
  // if (x.length() != m || g.length() != m) {
  //   cout << "Verification not passed: bad length\n";
  //   return false;
  // }
  //
  // ZZ_p::init(q);
  //
  // // test of digest e
  // ZZ_p hashzz;
  // hash_ZZp(hashzz,g,x,dleq->a);
  // if (dleq->e != hashzz) {
  //   cout << "Verification not passed: bad digest\n";
  //   return false;
  // }
  //
  // // test of values of a_i = g_i^z * x_i^e
  // ZZ_pPush push(p);
  // ZZ_p tmp1, tmp2, tmp3;
  // for (int i = 0; i < m; i++) {
  //   power(tmp2, g[i], rep(dleq->z));
  //   power(tmp3, x[i], rep(dleq->e));
  //   mul(tmp1, tmp2, tmp3);
  //   if (dleq->a[i] != tmp1) {
  //     cout << "Verification not passed: bad a_i\n";
  //     return false;
  //   }
  // }
  //
  // cout << "Verification passed\n";
  // return true;
// }

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////// LDEI //////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// in ZZq
// proof_t* ldei_prove(const ZZ& q, const ZZ& p, const Vec<ZZ_p>& g,
// const Vec<ZZ_p>& alpha, const long k, const Vec<ZZ_p>& x, const ZZ_pX& P) {
//   ZZ_p::init(q);
//   // verify entries
//   int m = g.length();
//   if (alpha.length() != m || x.length() != m || k > m)
//     return NULL;
//
//
//   //create the vector a
//   ZZ_pX R = random_ZZ_pX(k+1);
//   Vec<ZZ_p> R_eval = eval(R, alpha);
//
//   ZZ_p::init(p);
//
//   // create proof structure
//   proof_t *ld = proof_alloc(m);
//
//   for (int i = 0; i < m; i++) {
//     power(ld->a[i], g[i], rep(R_eval[i]));
//   }
//
//   // in ZZq
//   ZZ_pPush push(q);
//
//   // create e the digest of the hash function
//   hash_ZZp(ld->e, x, ld->a);
//
//   // create z the polynomial of ld
//   ZZ_pX tmp = ld->e * P;
//   ld->z = tmp + R;
//
//   return ld;
// }

// bool ldei_verify(const proof_t* ld, const ZZ& q, const ZZ& p,
// const Vec<ZZ_p>& g, const Vec<ZZ_p>& alpha, const long k, const Vec<ZZ_p>& x){
  // // test of entries
  // if (!ld)
  //   return false;
  // int m = ld->a.length();
  // if (alpha.length() != m || x.length() != m || g.length() != m) {
  //   cout << "Verification not passed: bad length\n";
  //   return false;
  // }
  //
  // // test degree polynomial
  // if (deg(ld->z) > k) {
  //   cout << "Verification not passed: bad degree of z\n";
  //   return false;
  // }
  //
  // ZZ_p::init(q);
  //
  // // test of e = hash(x,a)
  // ZZ_p hashzz;
  // hash_ZZp(hashzz, x, ld->a);
  // if (ld->e != hashzz) {
  //   cout << "Verification not passed: bad digest\n";
  //   return false;
  // }


  // test of x_i^e * a_i = g_i^z(alpha_i)
//   Vec<ZZ_p> zi = eval(ld->z, alpha);
//   ZZ_pPush push(p);
//   ZZ_p tmp1, tmp2, tmp3;
//   for (int i = 0; i < m; i++) {
//     power(tmp2, g[i], rep(zi[i]));
//     power(tmp3, x[i], rep(ld->e));
//     mul(tmp1, tmp3, ld->a[i]);
//     if (tmp2 != tmp1) {
//       cout << "Verification not passed: bad a_i\n";
//       return false;
//     }
//   }
//
//   cout << "Verification passed\n";
//   return true;
// }


////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// local LDEI ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// in ZZq
bool localldei(const ZZ& q, const ZZ& p, const Vec<ZZ_p>& alpha, const long k, const Vec<ZZ_p>& x,
const long m) {
  ZZ_p::init(q);
  // computation of vector u
  Vec<ZZ_p> u;
  u.SetLength(m);
  ZZ_p prod, tmp;
  for (int i = 0; i < m; i++) {
    prod = 1;
    for (int l = 0; l < m; l++)
      if (l != i) {
        sub(tmp, alpha[i], alpha[l]);
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
    eval(tmp, P, alpha[i]);
    mul(v[i], u[i], tmp);
  }
  // verification
  ZZ_pPush push(p);
  prod = 1;
  for (int i = 0; i < m; i++) {
    power(tmp, x[i], rep(v[i]));
    prod *= tmp;
  }
  v.kill();
  return (IsOne(prod));
}
