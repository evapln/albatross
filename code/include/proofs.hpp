#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/vector.h>
#include <NTL/ZZ_pX.h>

using namespace std;
using namespace NTL;


class LDEI {
  private:
    Vec<ZZ_p> a;
    ZZ_p e;
    ZZ_pX z;
  public:
    LDEI();
    ~LDEI() {}
    void print();
    void prove(const ZZ& q, const ZZ& p, const Vec<ZZ_p>& g, const Vec<ZZ_p>& alpha,
      const long k, const Vec<ZZ_p>& x, const ZZ_pX& P);
    bool verify(const ZZ& q, const ZZ& p, const Vec<ZZ_p>& g, const Vec<ZZ_p>& alpha,
      const long k, const Vec<ZZ_p>& x);
};


class DLEQ {
  private:
    Vec<ZZ_p> a;
    ZZ_p e;
    ZZ z;
  public:
    DLEQ();
    ~DLEQ() {}
    void print();
    void prove(const ZZ& q, const ZZ& p, const Vec<ZZ_p>& g, const Vec<ZZ_p>& x,
      ZZ_p& alpha);
    bool verify(const ZZ& q, const ZZ& p, const Vec<ZZ_p>& g, const Vec<ZZ_p>& x);
};


////////////////////////////////////////////////////////////////////////////////
///////////////////////////// proof structure //////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/* structure proof: (a_1,...,a_m,e,z)
- q prime, order of the group G
- g vector of m generators of G
- alpha vector of m pairwise distinct element of Z_q
- 1 <= k < m
- x vector of m elemenet of G
*/

// typedef struct proof_t proof_t;
//
// proof_t *proof_alloc(const int m);
//
// /* v only for test v */
// proof_t* proof_copy(proof_t* src, int m);
// void proof_set_a(proof_t* ld, Vec<ZZ_p>& a);
// void proof_set_e(proof_t* ld, ZZ_p& e);
// void proof_set_z(proof_t* ld, ZZ_pX& z);
// /* ^ only for test ^ */
//
// void proof_free(proof_t *ld);
//
// void proof_print(proof_t *ld);
//
// ////////////////////////////////////////////////////////////////////////////////
// //////////////////////////////////// DLEQ //////////////////////////////////////
// ////////////////////////////////////////////////////////////////////////////////
//
// proof_t* dleq_prove(const ZZ& q, const ZZ& p, const Vec<ZZ_p>& g,
// const Vec<ZZ_p>& x, ZZ_p& alpha);
//
// bool dleq_verify(const proof_t* dleq, const ZZ& q, const ZZ& p,
// const Vec<ZZ_p>& g, const Vec<ZZ_p>& x);
//
// ////////////////////////////////////////////////////////////////////////////////
// //////////////////////////////////// LDEI //////////////////////////////////////
// ////////////////////////////////////////////////////////////////////////////////
//
// proof_t* ldei_prove(const ZZ& q, const ZZ& p, const Vec<ZZ_p>& g, const Vec<ZZ_p>& alpha,
//   const long k, const Vec<ZZ_p>& x, const ZZ_pX& P);
//
// bool ldei_verify(const proof_t* ld, const ZZ& q, const ZZ& p, const Vec<ZZ_p>& g,
//   const Vec<ZZ_p>& alpha, const long k, const Vec<ZZ_p>& x);

////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// local LDEI ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
bool localldei(const ZZ& q, const ZZ& p, const Vec<ZZ_p>& alpha, const long k, const Vec<ZZ_p>& x,
  const long m);


/* time for one verification with q and vectors of size 1024 -> about 1.5s */
