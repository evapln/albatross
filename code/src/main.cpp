#include "ffte.hpp"
#include "pvss.hpp"
#include "proofs.hpp"
#include "hash.hpp"


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
  //////////////////////////////////////////////////////////////////// test ffte
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

  //////////////////////////////////////////////////////////////////// test pvss
  // cout << "pvss" << endl;
  // pvss_test();


  //////////////////////////////////////////////////////////////////// test DLEQ
  // cout << "DLEQ" << endl;
  // ZZ q = GenGermainPrime_ZZ(10);
  // ZZ p = 2 * q + 1;
  // int len = 5;
  // cout << "q = " << q << "  p = " << p << endl;
  // ZZ_p::init(q);
  // ZZ_p alpha = random_ZZ_p();
  // cout << "alpha = " << alpha << endl;
  // ZZ_p::init(p);
  // Vec<ZZ_p> g, x;
  // g.SetLength(len);
  // x.SetLength(len);
  // for (int i = 0; i < len; i++) {
  //   random(g[i]);
  //   power(x[i],g[i],rep(alpha));
  // }
  // cout << "x = " << x << endl << "g = " << g << endl;
  // DLEQ dl;
  // dl.prove(q,p,g,x,alpha);
  // dl.print();
  // cout << dl.verify(q,p,g,x) << endl;


  ////////////////////////////////////////////////////////////// test Local LDEI
  // cout << "Local LDEI vec" << endl;
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
  // verif1 = localldei(q, p, a, k, x, m);
  // time1 = clock() - rec;
  //
  // rec = clock();
  // verif2 = localldei(q, p, a, k, xf, m);
  // time2 = clock() - rec;
  //
  // cout << "It's " << verif1 << " in " << (float)time1/CLOCKS_PER_SEC << " second(s).\n";
  // cout << "It's " << verif2 << " in " << (float)time2/CLOCKS_PER_SEC << " second(s).\n";

  // cout << "Local LDEI int*" << endl;
  // long m = 1024, k = 764;
  // ZZ q = GenGermainPrime_ZZ(1024);
  // ZZ p = 2 * q + 1;
  // cout << "q = " << q << "  p = " << p << endl;
  // ZZ_p::init(q);
  // ZZ_pX P;
  // random(P,k+1);
  // int a[m];
  // Vec<ZZ_p> pi;
  // pi.SetLength(m);
  // prng_init(time(NULL));
  // for (int i = 0; i < m; i ++) {
  //   a[i] = rand() % 100;
  //   pi[i] = eval(P,ZZ_p(a[i]));
  // }
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
  // verif1 = localldei(q, p, a, k, x, m);
  // time1 = clock() - rec;
  //
  // rec = clock();
  // verif2 = localldei(q, p, a, k, xf, m);
  // time2 = clock() - rec;
  //
  // cout << "It's " << verif1 << " in " << (float)time1/CLOCKS_PER_SEC << " second(s).\n";
  // cout << "It's " << verif2 << " in " << (float)time2/CLOCKS_PER_SEC << " second(s).\n";

  //////////////////////////////////////////////////////////////////// test LDEI
  // cout << "LDEI" << endl;
  // long m = 1024, k = 9;
  // ZZ q = GenGermainPrime_ZZ(1024);
  // ZZ p = 2 * q + 1;
  // ZZ_p::init(q);
  // Vec<ZZ_p> alpha = vect_rand_dist(m);
  // ZZ_pX P = random_ZZ_pX(k+1);
  // Vec<ZZ_p> P_eval = eval(P,alpha);
  //
  // ZZ_p::init(p);
  // Vec<ZZ_p> gens;
  // generators(gens, q, m);
  // if (gens.length() != m)
  //   return EXIT_FAILURE;
  // Vec<ZZ_p> h;
  // h.SetLength(m);
  // for (int i = 0; i < m; i++)
  //   power(h[i], gens[i], 2);
  // Vec<ZZ_p> x;
  // x.SetLength(m);
  // for (int i = 0; i < m; i++)
  //   power(x[i],h[i],rep(P_eval[i]));
  //
  // clock_t rec;
  // rec = clock();
  // LDEI ld;
  // ld.prove(q, p, h, alpha, k, x, P);
  // cout << "should be true : ";
  // rec = clock() - rec;
  // ld.verify(q, p, h, alpha, k, x);
  // cout << (float)rec/CLOCKS_PER_SEC << "s\n";


  //////////////////////////////////////////////////////////////////// test hash
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

  //////////////////////////////////////////////////////////////// test mat bool
  // int n = 128, k = 20;
  // Mat<bool> M;
  // M.SetDims(k,n);
  // prng_init(time(NULL));
  // for (int i = 0; i < k; i++) {
  //   for (int j = 0; j < n; j++)
  //     M[i][j] = rand() % 2;
  // }
  //
  // for (int i = 0; i < k; i++) {
  //   for (int j = 0; j < n; j++)
  //     cout << M[i][j] << " ";
  //   cout << endl;
  // }
  //
  ///////////////////////////////////////////////////////// test mat mat_set_row
  // Vec<bool> row;
  // row.SetLength(n);
  // for (int i = 0; i < n; i++)
  //   row[i] = rand() % 2;
  // cout << "we set this row in pos 4 : " << row << endl;
  // mat_set_row(M, 4, row);
  //
  // for (int i = 0; i < k; i++) {
  //   for (int j = 0; j < n; j++)
  //     cout << M[i][j] << " ";
  //   cout << endl;
  // }


  ////////////////////////////////////////////////////////////// tes vec_is_zero
  // bool test1, test2;
  // test1 = vec_is_zero(row);
  // Vec<bool> zero;
  // zero.SetLength(n);
  // for (int i = 0; i < n; i++)
  //   zero[i] = 0;
  // test2 = vec_is_zero(zero);
  // cout << "should be 0 : " << test1 << endl;
  // cout << "should be 1 : " << test2 << endl;
  //
  // test vec_is_in_mat
  // Vec<bool> r0, r;
  // r0.SetLength(n);
  // r.SetLength(n);
  // for (int i = 0; i <n; i++) {
  //   r0[i] = M[2][i];
  //   r[i] = rand() %2 ;
  // }
  // for (int i = 0; i < k; i++) {
  //   for (int j = 0; j < n; j++)
  //     cout << M[i][j] << " ";
  //   cout << endl;
  // }
  // cout << "r0 " << r0 << " ";
  // cout << vec_is_in_mat(M,r0) << endl;
  // cout << "r " << r << " ";
  // cout << vec_is_in_mat(M,r) << endl;
  //
  // // test mat gen and space_add_row
  // Mat<bool> mat;
  // cout << mat_gen(mat,n,k) << endl;
  // for (int i = 0; i < k; i++) {
  //   for (int j = 0; j < n; j++)
  //     cout << mat[i][j] << " ";
  //   cout << endl;
  // }


  ///////////////////////////////////////////////////////////////// test mul_mat
  // ZZ q = GenGermainPrime_ZZ(n);
  // ZZ_p::init(q);
  // Vec<ZZ_p> hhat;
  // Vec<ZZ_p> h;
  // h.SetLength(n);
  // for (int i = 0; i < n; i++)
  //  random(h[i]);
  // cout << "q = " << q << "\nh " << h << endl;
  // mul_mat(hhat,mat,n,k,h);
  // if (hhat.length() == 0)
  //   return EXIT_FAILURE;
  // cout << "hhat " << hhat << endl;

  /////////////////////////////////////////////////// comparision FFTE / mul_mat
  int n = 2048, k = 16, t = 10, l = n-2*t;
  ZZ p,q;
  findprime(k,l,q,p);
  ZZ w;
  rootunity(w,n,q);
  // ZZ q = GenGermainPrime_ZZ(128);
  // cout << q << endl;
  ZZ_p::init(p);
  Mat<ZZ_p> M;
  clock_t rec, tsame = 0, tdif = 0, tffte = 0;
  int m = 10;
  Vec<ZZ_p> vech;
  Vec<ZZ_p> hhat;
  Vec<ZZ_p> f;
  vech.SetLength(n);
  ZZ_p h {4};


  Vec<ZZ_p> L, coef;
  L.SetLength(n);
  coef.SetLength(n);

  // test time m computations with the same matrix
  mat_gen(M,n,k);
  for (int i = 0; i < m; i++) {
    // new vector
    for (int j = 0; j < n; j++)
     random(vech[i]);

    rec = clock();
    mul_mat(hhat,M,n,k,vech);
    tsame += clock() - rec;
  }

  // test time m computations with the same matrix
  for (int i = 0; i < m; i++) {
    // new matrix
    mat_gen(M,n,k);
    // new vector
    for (int j = 0; j < n; j++)
     random(vech[i]);

    rec = clock();
    mul_mat(hhat,M,n,k,vech);
    tdif += clock() - rec;
  }

  // test ffte
  for (int i = 0; i < m; i++) {
    // new vector
    for (int i = 0; i < n; i++) {
      random(coef[i]);
      power(L[i],h,rep(coef[i]));
    }

    rec = clock();
    FFTE(f, n, L, w, q);
    tffte += clock() - rec;
  }

  cout << "same matrix, " << m << " computations : " << (float)tsame/CLOCKS_PER_SEC << " second(s).\n";
  cout << "different matrices, " << m << " computations : " << (float)tdif/CLOCKS_PER_SEC << " second(s).\n";
  cout << "ffte, " << m << " computations : " << (float)tffte/CLOCKS_PER_SEC << " second(s).\n";




  /////////////////////////////////////////////////////////////////////// return
  return EXIT_SUCCESS;
}
