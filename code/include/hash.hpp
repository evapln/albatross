#include <NTL/ZZ_p.h>

using namespace std;
using namespace NTL;

/* convert a string in a ZZ_p */
void string_to_ZZp(ZZ_p& output, const string& input);

/* convert a ZZ_p in a string */
string ZZp_to_string(const Vec<ZZ_p>& x, const Vec<ZZ_p>& a);

/* hash a message with sha3_512 and convert it in a ZZ_p */
string sha3_512(const string msg);

/* hashzz = hash(x,a) */
void hash_ZZp(ZZ_p& hashzz, const Vec<ZZ_p>& x, const Vec<ZZ_p>& a);
