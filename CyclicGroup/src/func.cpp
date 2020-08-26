#include <cstdlib>
#include <stdbool.h>

#include "func.hpp"


void prng_init(const unsigned int seed) {
  static bool seed_init = false;
  if (! seed_init) {
    srand(seed);
    seed_init = true;
  }
}

void findprime(ZZ& q, ZZ& p, const long k, const long l) {
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

void generator(ZZ_p& g, const ZZ& q) {
  ZZ_p po;
  for (int i = 2; i < 2*q+1; i++) {
    power(po,ZZ_p(i),2);
    if (po == 1)
      continue;
    power(po,ZZ_p(i),q);
    if (IsOne(po))
      continue;
    g = ZZ_p(i);
    return;
  }
  return;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////// gestion of vectors and matricse /////////////////////////
////////////////////////////////////////////////////////////////////////////////

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
  // for all columns
  for (int j = 0; j < col; ++j) {
    // for all elements in the column
    for (int i = j; i < row; ++i) {
      // put an non zero element on top of each column
      if (i != j && !IsZero(M[i][j])) {
        mat_exchange_row(M, j, i);
        break;
      }
    }
    // set  the rest of the column at zero
    for (int i = j + 1; i < row; ++i) {
      if (!IsZero(M[i][j])) {
        inv(inve, M[j][j]);
        mul(coef, M[i][j], inve);
        if (!IsZero(coef))
          mat_add_row(M, i, j, -coef);
      }
    }
  }
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
