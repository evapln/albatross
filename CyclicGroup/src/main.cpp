/* Implementation of ppvss over cyclic group */

#include "randomness_extraction.hpp"
#include "pvss.hpp"
#include "proofs.hpp"
#include "hash.hpp"
#include "func.hpp"

#include <cstdlib>
#include <ctime>
#include <iostream>
#include <getopt.h>
#include <NTL/ZZ_pX.h>

int main(int argc, char *argv[]) {
  int optc;
  int n = 1024,size = 1024;
  bool ppvss = false, comp = false, ffte = false;

  /* long options struct */
  const struct option long_opts[] =
  {
    {"ppvss",                  optional_argument, NULL, 'p'},
    {"comparison",             no_argument,       NULL, 'c'},
    {"ffte",                   optional_argument, NULL, 'f'},
    {"number_of_participants", required_argument, NULL, 'n'},
    {"help",                   no_argument,       NULL, 'h'},
    {NULL,                     0,                 NULL,  0 },
  };

  const char *opt = "p::cf::n:h";
  while ((optc = getopt_long(argc, argv, opt, long_opts, NULL)) != -1) {
    switch (optc) {
      case 'p':
        if(optarg != NULL)
          size = atoi(optarg);
        ppvss = true;
        break;
      case 'c':
        comp = true;
        break;
      case 'f':
        if(optarg != NULL)
          size = atoi(optarg);
        ffte = true;
        break;
      case 'n':
        n = atoi(optarg);
        break;
      case 'h':
        cout << "Usage :\tALBATROSS -p [SIZE_OF_Q|-n NUMBER_OF_PARTICIPANTS|-h]\n"
                "\tALBATROSS -c[-n NUMBER_OF_PARTICIPANTS|-h]\n"
                "\tALBATROSS -f [SIZE_OF_Q|-n NUMBER_OF_PARTICIPANTS|-h]\n\n"
                "ALBATROSS: publicly AttestabLe BATched Randomness based On Secret Sharing \n\n"
                "-p, --ppvss\t\t\trun the ppvss scheme with q of size the argument or 1024 bits if not specified\n"
                "-c, --comparison\t\tcompare the two method of extraction\n"
                "-f, --ffte\t\t\trun the ffte with q of size the argument or 1024 bits if not specified\n"
                "-n, --number_of_participants\tset n, without this option, n = 1024\n"
                "-h, --help\t\t\tdisplay this help\n\n";
        return EXIT_SUCCESS;
      default:
        cout << "error: invalid option!\n";
        return EXIT_FAILURE;
    }
  }

  if ((ppvss && comp) || (ffte && comp) || (ffte && ppvss)){
    cout << "Please, enter one mode at a time\n";
    return EXIT_FAILURE;
  }

  if (!ppvss && !comp && !ffte){
    cout << "Please, enter a mode (p, r, f or h)\n";
    return EXIT_FAILURE;
  }

  if (ppvss) {
    ///////////////////////////////////////////////////////////////// test pvss
    cout << "execution of ppvss\n\n";
    pvss_test(n,size);
  }

  if (comp) {
    /////////////////////////////////////////////// comparision FFTE / mul_mat
    // int n = 2048;
    cout << "comparison between ffte and 0-1 matrix\n\n";
    cout << n << endl;
    int k = n/128;
    ZZ p,q;
    findprime(q,p,100,28);
    ZZ w;
    rootunity(w,n,q);
    ZZ_p::init(p);
    Mat<ZZ_p> M;
    clock_t rec, tsame = 0, tdif = 0, tffte = 0;
    int m = 1000;
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
  }

  if (ffte) {
    cout << "execution of ffte\n\n";
    long k = 128;
    long l = size - k;
    ZZ p,q;
    findprime(q,p,k,l);
    cout << "q = " << q << "  p = " << p << endl;
    ZZ w;
    rootunity(w,n,q);
    cout << "w = " << w << endl;
    ZZ_p::init(p);
    ZZ_p h {4};
    Vec<ZZ_p> L, coef;
    L.SetLength(n);
    coef.SetLength(n);
    for (int i = 0; i < n; i++) {
      random(coef[i]);
      power(L[i],h,rep(coef[i]));
    }
    Vec<ZZ_p> f;
    clock_t FFTE_time = clock();
    FFTE(f,n,L,w,q);
    FFTE_time = clock() - FFTE_time;
    cout << "f = " << f << endl;
    cout << "correct ? " << test(f,coef,h,w,n,q) << endl;
    cout << "time: " << (float)FFTE_time/CLOCKS_PER_SEC << "s" << endl;
  }

  /////////////////////////////////////////////////////////////////////// return
  return EXIT_SUCCESS;
}
