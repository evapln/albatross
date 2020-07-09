#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/vector.h>

using namespace std;
using namespace NTL;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// LDEI ////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/* structure ldei: (a_1,...,a_m,e,z)
- q prime, order of the group G
- g vector of m generators of G
- alpha vector of m pairwise distinct element of Z_q
- 1 <= k < m
- x vector of m elemenet of G
*/

typedef struct ldei_t ldei_t;

ldei_t *ldei_alloc(const int m);


ldei_t* ldei_copy(ldei_t* src, int m);
void ldei_set_a(ldei_t* ld, Vec<ZZ_p>& a);
void ldei_set_e(ldei_t* ld, ZZ_p& e);
void ldei_set_z(ldei_t* ld, ZZ_pX& z);

void ldei_free(ldei_t *ld);

void ldei_print(ldei_t *ld);

ldei_t* ldei_prove(const ZZ& q, const ZZ& p, const Vec<ZZ_p>& g, const Vec<ZZ_p>& alpha,
  const long k, const Vec<ZZ_p>& x, const ZZ_pX& P);

bool ldei_verify(const ldei_t* ld, const ZZ& q, const ZZ& p, const Vec<ZZ_p>& g,
  const Vec<ZZ_p>& alpha, const long k, const Vec<ZZ_p>& x);

////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// local LDEI ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

bool localldei(const Vec<ZZ_p>& a, const long k, const Vec<ZZ_p>& x,
  const long m, const ZZ& p);



/* time for one verification with q and vectors of size 1024 -> about 1.5s */
