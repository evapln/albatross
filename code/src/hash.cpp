#include <cryptopp/sha3.h>
#include <sstream>
#include <string>
#include "../include/hash.hpp"

using namespace CryptoPP;

void string_to_ZZp(ZZ_p& output, const string& input) {
  size_t len = input.length();
  output = ZZ_p(0);
  for (size_t i = 0; i < len; i++) {
    const unsigned char c = input[i];
		mul(output,output, 1000);
		add(output, output, int(c));
  }
}

string ZZp_to_string(const Vec<ZZ_p>& x, const Vec<ZZ_p>& a) {
  stringstream vecx, veca;
  string vecxx, vecaa, s = "";
  char c;
  vecx << x;
  vecxx = vecx.str();
  for (size_t i = 0; i < vecxx.length(); i++) {
    c = vecxx[i];
    switch(c) {
      case '[':
      case ']':
      case ' ':
        break;
      default:
        s.push_back(c);
    }
  }
  veca << a;
  vecaa = veca.str();
  for (size_t i = 0; i < vecaa.length(); i++) {
      c = vecaa[i];
      switch(c) {
        case '[':
        case ']':
        case ' ':
          break;
        default:
          s.push_back(c);
      }
    }
  return  s;
}


string sha3_512(const string msg) {
	string digest;
	SHA3_512 hash;
	hash.Update((const byte*)msg.data(), msg.size());
	digest.resize(hash.DigestSize());
	hash.Final((byte*)&digest[0]);
	return digest;
}

void hash_ZZp(ZZ_p& hashzz, const Vec<ZZ_p>& x, const Vec<ZZ_p>& a) {
  string msg = ZZp_to_string(x,a);
  cout << "msg = " << msg << endl;
	string digest = sha3_512(msg);
	string_to_ZZp(hashzz, digest);
}
