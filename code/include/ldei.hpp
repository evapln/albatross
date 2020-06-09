#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/vector.h>

using namespace std;
using namespace NTL;

bool localldei(const Vec<ZZ_p>& a, const long k, const Vec<ZZ_p>& x,
  const long m, const ZZ& p);

/* time for one verification with q and vectors of size 1024 -> about 1.5s */
