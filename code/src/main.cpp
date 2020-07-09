#include "../include/ffte.hpp"
#include "../include/pvss.hpp"
#include "../include/ldei.hpp"
#include "../include/hash.hpp"


#include <iostream>
#include <cstdlib>
#include <ctime>
#include <NTL/ZZ_pX.h>

bool IsIn(const Vec<ZZ_p>& v, const int n, const ZZ_p val) {
    for (int i = 0; i < n; i++)
      if (v[i] == val)
        return true;
    return false;
}

Vec<ZZ_p> vect_rand_dist(const int n) {
  Vec<ZZ_p> v;
  v.SetLength(n);
  ZZ_p tmp;
  int len = 0;
  while (len < n) {
    random(tmp);
    if (!IsIn(v, len, tmp)) {
      v[len] = tmp;
      len++;
    }
  }
  return v;
}

// in ZZp
void generators(Vec<ZZ_p>& g, const ZZ& q, const int m) {
  g.SetLength(m);
  ZZ_p po;
  int ind = 0;
  for (int i = 2; i < 2*q+1; i++) {
    power(po,ZZ_p(i),2);
    if (IsOne(po))
      continue;
    power(po,ZZ_p(i),q);
    if (IsOne(po))
      continue;
    g[ind] = ZZ_p(i);
    ind++;
    if (ind >= m)
      return;
  }
  if (ind < m)
    g.kill();
  return;
}

string string_to_hex(const string& input) {
  static const char* const lut = "0123456789ABCDEF";
  size_t len = input.length();
  string output;
  output.reserve(2 * len);
  for (size_t i = 0; i < len; ++i) {
    const unsigned char c = input[i];
    output.push_back(lut[c >> 4]);
    output.push_back(lut[c & 15]);
  }
	//cout<<" output"<<output<<endl;
  return output;
}

int main(void) {
  // // test ffte
  // cout << "ffte" << endl;
  // long n = 1024;
  // long k = 102;
  // long l = 1024 - k;
  // ZZ p,q;
  // findprime(k,l,q,p);
  // cout << "q = " << q << "  p = " << p << endl;
  // ZZ w;
  // rootunity(w,n,q);
  // cout << "w = " << w << endl;
  // ZZ_p::init(p);
  // ZZ_p h {4};
  // Vec<ZZ_p> L, coef;
  // L.SetLength(n);
  // coef.SetLength(n);
  // for (int i = 0; i < n; i++) {
  //   random(coef[i]);
  //   power(L[i],h,rep(coef[i]));
  // }
  // Vec<ZZ_p> f;
  // clock_t FFTE_time = clock();
  // FFTE(f,n,L,w,q);
  // FFTE_time = clock() - FFTE_time;
  // cout << "f = " << f << endl;
  // cout << "correct ? " << test(f,coef,h,w,n,q) << endl;
  // cout << "time: " << (float)FFTE_time/CLOCKS_PER_SEC << "s" << endl;

  // // test pvss
  // cout << "pvss" << endl;
  // pvss();

  // test Local LDEI
  // cout << "Local LDEI" << endl;
  // long m = 1024, k = 764;
  // ZZ q = GenGermainPrime_ZZ(1024);
  // ZZ p = 2 * q + 1;
  // cout << "q = " << q << "  p = " << p << endl;
  // ZZ_p::init(q);
  // Vec<ZZ_p> a = vect_rand_dist(m);
  // ZZ_pX P;
  // random(P,k+1);
  // Vec<ZZ_p> pi = eval(P,a);
  //
  // ZZ_p::init(p);
  // ZZ_p g;
  // generator(g,p);
  // ZZ_p h = power(g,2);
  // cout << g << endl << h <<  endl;
  // Vec<ZZ_p> x;
  // x.SetLength(m);
  // for (int i = 0; i < m; i++)
  //   power(x[i],h,rep(pi[i]));
  //
  // Vec<ZZ_p> xf;
  // xf.SetLength(m);
  // for (int i = 0; i < m; i++)
  //   random(xf[i]);
  //
  // ZZ_p::init(q);
  // clock_t rec, time1, time2;
  // bool verif1, verif2;
  //
  // rec = clock();
  // verif1 = localldei(a, k, x, m, p);
  // time1 = clock() - rec;
  //
  // rec = clock();
  // verif2 = localldei(a, k, xf, m, p);
  // time2 = clock() - rec;
  //
  // cout << "It's " << verif1 << " in " << (float)time1/CLOCKS_PER_SEC << " second(s).\n";
  // cout << "It's " << verif2 << " in " << (float)time2/CLOCKS_PER_SEC << " second(s).\n";

  // test LDEI
  cout << "LDEI" << endl;
  long m = 1024, k = 9;
  ZZ q = GenGermainPrime_ZZ(1024);
  ZZ p = 2 * q + 1;
  ZZ_p::init(q);
  Vec<ZZ_p> alpha = vect_rand_dist(m);
  ZZ_pX P = random_ZZ_pX(k+1);
  Vec<ZZ_p> P_eval = eval(P,alpha);

  ZZ_p::init(p);
  Vec<ZZ_p> gens;
  generators(gens, q, m);
  if (gens.length() != m)
    return EXIT_FAILURE;
  Vec<ZZ_p> h;
  h.SetLength(m);
  for (int i = 0; i < m; i++)
    power(h[i], gens[i], 2);
  Vec<ZZ_p> x;
  x.SetLength(m);
  for (int i = 0; i < m; i++)
    power(x[i],h[i],rep(P_eval[i]));

  clock_t rec;
  rec = clock();
  ldei_t *ld = ldei_prove(q, p, h, alpha, k, x, P);
  rec = clock() - rec;
  cout << (float)rec/CLOCKS_PER_SEC << "s\n";
  cout << "should be true : ";
  ldei_verify(ld, q, p, h, alpha, k, x);

  ldei_t* copy1 = ldei_copy(ld, m);
  ldei_t* copy2 = ldei_copy(ld, m);
  ldei_t* copy3 = ldei_copy(ld, m);
  ldei_t* copy4 = ldei_copy(ld, m);
  ldei_free(ld);
  ZZ_p::init(q);


  ZZ_p ef = random_ZZ_p();
  ldei_set_e(copy1, ef);
  cout << "should be false from hash: ";
  rec = clock();
  ldei_verify(copy1, q, p, h, alpha, k, x);
  rec = clock() - rec;
  cout << (float)rec/CLOCKS_PER_SEC << "s\n";
  ldei_free(copy1);

  ZZ_pX zf = random_ZZ_pX(k+1);
  ldei_set_z(copy2, zf);
  cout << "should be false from ai: ";
  rec = clock();
  ldei_verify(copy2, q, p, h, alpha, k, x);
  rec = clock() - rec;
  cout << (float)rec/CLOCKS_PER_SEC << "s\n";
  ldei_free(copy2);

  Vec<ZZ_p> af;
  af.SetLength(m);
  for (int i = 0 ; i < m; i++)
    random(af[i]);
  ldei_set_a(copy3, af);
  cout << "should be false from hash: ";
  rec = clock();
  ldei_verify(copy3, q, p, h, alpha, k, x);
  rec = clock() - rec;
  cout << (float)rec/CLOCKS_PER_SEC << "s\n";
  ldei_free(copy3);

  af.SetLength(m+1);
  for (int i = 0 ; i < m+1; i++)
    random(af[i]);
  ldei_set_a(copy4, af);
  cout << "should be false from lengths: ";
  rec = clock();
  ldei_verify(copy4, q, p, h, alpha, k, x);
  rec = clock() - rec;
  cout << (float)rec/CLOCKS_PER_SEC << "s\n";
  ldei_free(copy4);




  // // test file hash
  // cout << "hash" << endl;
  // long k = 102;
  // long l = 1024 - k;
  // ZZ p,q;
  // findprime(k,l,q,p);
  // cout << q << endl;
  // ZZ_p::init(q);
  // ZZ_p zz_s2z, out;
  // Vec<ZZ_p> x, a;
  // x.SetLength(5);
  // for (int i = 0; i < 5; i++) x[i] = ZZ_p(i+11);
  // a.SetLength(4);
  // for (int i = 0; i < 4; i++) a[i] = ZZ_p(i+16);
  // string s_s2z, s_z2s, in_sha, out_sha;
  // // s2z
  // cout << "string_to_ZZp" << endl;
  // s_s2z = "test";
  // string_to_ZZp(zz_s2z, s_s2z);
  // cout << "sould be " << s_s2z << " : " << zz_s2z << endl;
  // // z2s
  // cout << "ZZP_to_string" << endl;
  // s_z2s = ZZp_to_string(x, a);
  // cout << "should be " << x << a << " : " << s_z2s << endl;
  // // sha
  // cout << "sha3_512" << endl;
  // in_sha = "111213141516171819";
  // out_sha = sha3_512(in_sha);
  // cout << "hash of " << in_sha << " is " << string_to_hex(out_sha) << endl;
  // // hash_ZZp
  // cout << "hash_ZZp" << endl;
  // hash_ZZp(out, x, a);
  // cout << out << endl;

  // return
  return EXIT_SUCCESS;
}
