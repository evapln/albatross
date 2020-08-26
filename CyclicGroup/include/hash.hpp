#include <NTL/ZZ_p.h>

using namespace std;
using namespace NTL;

/* convert the input string in a ZZ_p */
void string_to_ZZp(ZZ_p& output, const string& input);

/* return the concatenation of [g],x and a as a string */
string ZZp_to_string(const Vec<ZZ_p>& x, const Vec<ZZ_p>& a);
string ZZp_to_string(const Vec<ZZ_p>& g, const Vec<ZZ_p>& x, const Vec<ZZ_p>& a);

/* hash a message with sha3_512 and return the digest */
string sha3_512(const string msg);

/* set hashzz to hash(concatenation of [g],x,a) */
void hash_ZZp(ZZ_p& hashzz, const Vec<ZZ_p>& x, const Vec<ZZ_p>& a);
void hash_ZZp(ZZ_p& hashzz, const Vec<ZZ_p>& g, const Vec<ZZ_p>& x, const Vec<ZZ_p>& a);
