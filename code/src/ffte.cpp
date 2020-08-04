#include <ctime>
#include <iostream>
#include <NTL/ZZ.h>
#include <math.h>
#include <stdbool.h>
#include <unistd.h>

#include "ffte.hpp"
#include "func.hpp"


using namespace std;

void findprime(const long k, const long l, ZZ& q, ZZ& p) {
  ZZ n, tmp;
  power2(n, k);
  long s = k % 2 - l % 2;
  power2(tmp, l); // tmp = 2**l
  q = (tmp + s)* n + 1; // q = n*(2**l+s) + 1
  p = 2 * q + 1;
  int r = 0;
  int bound = pow(10,8);
  while (r < bound) {
    if (ProbPrime(q) == 1 && ProbPrime(p) == 1)
      return;
    q += 3 * n;
    p += 6 * n;
    r++;
  }
}

// mod q
void rootunity(ZZ& w, const long n, const ZZ& q) {
  int i = 2;
  ZZ t = (q-1) / n;
  ZZ tmp;
  while (true) {
    PowerMod(w,ZZ(i),t,q);
    PowerMod(tmp, w, n/2,q);
    if (!IsOne(tmp))
      return;
    i++;
  }
}

// dans ZZp
void FFTE(Vec<ZZ_p>& f, const long n, const Vec<ZZ_p>& h, const ZZ& w, const ZZ& q) {
  f.SetLength(n);
  if (n == 1) {
    f[0] = h[0];
    return;
  }
  long m = n/2;
  Vec<ZZ_p> u; //vec de ZZp
  u.SetLength(m);
  Vec<ZZ_p> v; // vec de ZZp
  v.SetLength(m);
  ZZ_p hj, hjm, invhjm, tmp;
  ZZ wj, w2; // in ZZq
  for (int j = 0; j < m; j++) {
    hj = h[j];
    hjm = h[j+m];
    u[j] = hj*hjm;
    inv(invhjm,hjm);
    PowerMod(wj,w,j,q);
    mul (tmp, hj, invhjm);
    power(v[j], tmp, wj);
  }
  Vec<ZZ_p> uhat;
  uhat.SetLength(m);
  Vec<ZZ_p> vhat;
  vhat.SetLength(m);
  PowerMod(w2,w,2,q);
  FFTE(uhat,m,u,w2,q);
  FFTE(vhat,m,v,w2,q);
  for (int i = 0; i < m; i++) {
    f[i*2] = uhat[i];
    f[i*2+1] = vhat[i];
  }
}

// dans ZZq
void evalu(ZZ& sum, const Vec<ZZ_p>& L, const ZZ& x, const ZZ& q) {
  sum = 0;
  ZZ tmp, tmp1;
  for (int i = 0; i < L.length(); i++) {
    PowerMod(tmp, x, i, q);
    MulMod(tmp1, rep(L[i]), tmp, q);
    AddMod(sum, sum, tmp1, q);
  }
}

// dans ZZp
bool test(const Vec<ZZ_p>& W, const Vec<ZZ_p>& coef, const ZZ_p& h, const ZZ& w,
          const long n, const ZZ& q) {
  ZZ eval;
  ZZ wi;
  Vec<ZZ_p> U;
  U.SetLength(n);
  for (int i = 0; i < n; i++) {
    PowerMod(wi,w,i,q);
    evalu(eval, coef, wi, q);
    power(U[i],h, eval);
  }
  return (U == W);
}

int mat_set_row(Mat<ZZ_p>& M, const int i, const Vec<ZZ_p>& vec) {
  int n = M.NumCols();
  if (n != vec.length() || i >= M.NumRows())
    return 0;
  for (int j = 0; j < n; j++)
    M[i][j] = vec[j];
  return 1;
}

bool vec_is_zero(const Vec<ZZ_p>& vec) {
  for (int i = 0; i < vec.length(); i++)
    if (!IsZero(vec[i]))
      return false;
  return true;
}

// bool vec_is_in_mat(const Mat<bool>& M, const Vec<bool>& v) {
//   int n = M.NumCols();
//   if (n != v.length())
//     return false;
//   bool same;
//   for (int i = 0; i < M.NumRows(); i++) {
//     same = true;
//     for (int j = 0; j < n; j++)
//       if (v[j] != M[i][j]) {
//         same = false;
//         break;
//       }
//     if (same)
//       return true;
//   }
//   return false;
// }


// int space_add_row(Mat<bool>& S, const int i, const Vec<bool>& r) {
//   int n = S.NumCols(), dim = S.NumRows();
//   if (r.length() != n || pow(2,i+1) > dim) {
//     cout << "dim error\n" << dim << endl;
//     return 0;
//   }
//   Vec<bool> sum;
//   sum.SetLength(n);
//   for (int j = 0; j < dim; j++) {
//     for (int k = 0; k < n; k++)
//       sum[k] = r[k] ^ S[j][k];
//     mat_set_row(S,pow(2,i)+j, sum);
//   }
//   return 1;
// }


// int mat_gen_space(Mat<bool>& M, const int n, const int k) {
//   M.SetDims(k,n);
//   Mat<bool> space;
//   long dim = pow(2,k);
//   space.SetDims(dim,n);
//   for (int i = 0; i < dim; i++)
//     for (long j = 0; j < n; j++)
//       space[i][j] = 0;
//   Vec<bool> r;
//   r.SetLength(n);
//
//   // first row
//   // do {
//   //   for (int j = 0; j < n; j++)
//   //     r[j] = rand() % 2;
//   // } while (vec_is_zero(r));
//   // if (mat_set_row(M, 0, r) != 1 || space_add_row(space, 0, r) != 1)
//   //   return 0;
//
//   // next rows
//   bool in;
//   for (int i = 0; i < k; i++) {
//     in = true;
//     while (in) {
//       // new vector
//       for (int j = 0; j < n; j++)
//         r[j] = rand() % 2;
//       in = vec_is_in_mat(space,r);
//     }
//     if (mat_set_row(M,i,r) != 1 || space_add_row(space,i,r) != 1)
//       return 0;
//   }
//   // cout << "space\n";
//   // for (int i = 0; i < dim; i++)
//   //   cout << space[i] << endl;
//   return 1;
// }

void mat_exchange_row(Mat<ZZ_p>& M, const int row1, const int row2) {
  int row = M.NumRows();
  if (row1 > row || row2 > row)
    return;
  ZZ_p tmp;
  for (int i = 0; i < M.NumCols(); i++) {
    tmp = M[row1][i];
    M[row1][i] = M[row2][i];
    M[row2][i] = tmp;
  }
}

void mat_add_row(Mat<ZZ_p>& M, const int row1, const int row2, const ZZ_p& coef) {
  int row = M.NumRows();
  if (row1 > row || row2 > row)
    return;
  ZZ_p tmp;
  for (int i = 0; i < M.NumCols(); ++i) {
    mul(tmp, M[row2][i], coef);
    add(M[row1][i], M[row1][i], tmp);
  }
}

void mat_systematisation(Mat<ZZ_p>& M, const int row) {
  if (row > M.NumRows())
    return;
  int col = M.NumCols();
  ZZ_p coef, inve;
  // étape 1
  // pour toutes les colonnes
  for (int j = 0; j < col; ++j) {
    // parcours de chaque colonne
    for (int i = j; i < row; ++i) {
      // on met un coefficient non nul en haut de chaque colonne
      if (i != j && !IsZero(M[i][j])) {
        mat_exchange_row(M, j, i);
        break;
      }
    }
    // annulation du reste de la colonne
    for (int i = j + 1; i < row; ++i) {
      if (!IsZero(M[i][j])) {
        inv(inve, M[j][j]);
        mul(coef, M[i][j], inve);
        if (!IsZero(coef))
          mat_add_row(M, i, j, -coef);
      }

      // if (!IsZero(coef)) {
      //   matrix_add_row(M, i, j, -coef);
      // }
    }
  }

  // // on met tous les pivots à 1
  // // pour toutes les lignes
  // for (int i = 0; i < row; ++i) {
  //   // on parcours les lignes
  //   for (int j = 0; j < col; ++j) {
  //     // on trouve le premier coef non nul, le pivot
  //     if (matrix->mat[i][j] != 0) {
  //       matrix_mul_row(matrix, i, inv_Fq(matrix->mat[i][j]));
  //       break;
  //     }
  //   }
  // }
  // for (int i = row - 1; i >= 0; --i) {
  //   // parcours de la ligne
  //   int j = 0;
  //   while ((matrix->mat[i][j] == 0) && (j != col - 1))
  //     j += 1;
  //   for (int k = 0; k < i; ++k) {
  //     if (matrix->mat[k][j] != 0)
  //       matrix_add_row(matrix, k, i, -inv_Fq(matrix->mat[k][j]));
  //   }
  // }
}

bool is_ind(const Mat<ZZ_p>& M, const int row) {
  if (row > M.NumRows())
    return false;
  Mat<ZZ_p> copy = M;
  mat_systematisation(copy, row);
  for (int i = 0; i < row; i++)
    if(vec_is_zero(copy[i]))
      return false;
  return true;
}

Vec<ZZ_p> bin_rand_vec(const int n) {
  Vec<ZZ_p> vec;
  vec.SetLength(n);
  prng_init(time(NULL) + getpid());
  do {
    for (int i = 0; i < n; i++)
      vec[i] = ZZ_p(rand() % 2);
  } while (vec_is_zero(vec));
  return vec;
}

void bin_rand_vec(Vec<ZZ_p>& vec, const int n) {
  vec.SetLength(n);
  prng_init(time(NULL) + getpid());
  do {
    for (int i = 0; i < n; i++)
      vec[i] = ZZ_p(rand() % 2);
  } while (vec_is_zero(vec));
}

void mat_gen(Mat<ZZ_p>& M, const int n, const int k) {
  M.SetDims(k,n);
  Vec<ZZ_p> r;
  bool ind = false;
  while (!ind) {
    for (int i = 0; i < k; i++) {
      bin_rand_vec(r,n);
      mat_set_row(M,i,r);
    }
    ind = is_ind(M,k);
  }
}

void mul_mat(Vec<ZZ_p>& hhat, Mat<ZZ_p>& M, const int n, const int k, const Vec<ZZ_p>& h) {
  if (h.length() != n) {
    hhat.SetLength(0);
    return;
  }
  hhat.SetLength(k);
  for (int i = 0; i < k; i++) {
    hhat[i] = ZZ_p(1);
    for (int j = 0; j < n; j++) {
      if (M[i][j] == 1)
        hhat[i] *= h[j];
      }
  }
}
