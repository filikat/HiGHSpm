#include "CurtisReidScaling.h"

void product(const double* x, std::vector<double>& y,
             const std::vector<int>& ptr, const std::vector<int>& rows) {
  // Multiply by matrix E, i.e. matrix A with all entries equal to one
  // E * x = y
  int n = ptr.size() - 1;
  for (int col = 0; col < n; ++col) {
    for (int el = ptr[col]; el < ptr[col + 1]; ++el) {
      y[rows[el]] += x[col];
    }
  }
}

void product_transpose(const double* x, std::vector<double>& y,
                       const std::vector<int>& ptr,
                       const std::vector<int>& rows) {
  // Multiply by matrix E^T, i.e. matrix A^T with all entries equal to one
  // E^T * x = y
  int n = ptr.size() - 1;
  for (int col = 0; col < n; ++col) {
    for (int el = ptr[col]; el < ptr[col + 1]; ++el) {
      y[col] += x[rows[el]];
    }
  }
}

void mult_by_CR_matrix(const std::vector<double>& rhs, std::vector<double>& lhs,
                       const std::vector<double>& M,
                       const std::vector<double>& N,
                       const std::vector<int>& ptr,
                       const std::vector<int>& rows) {
  int n = N.size();
  int m = M.size();

  // split rhs
  const double* rho = &rhs.data()[0];
  const double* gamma = &rhs.data()[m];

  // compute E*gamma
  std::vector<double> Egamma(m);
  product(gamma, Egamma, ptr, rows);

  // compute E^T*rho
  std::vector<double> ETrho(n);
  product_transpose(rho, ETrho, ptr, rows);

  // populate lhs

  // 0...m-1
  for (int i = 0; i < m; ++i) lhs[i] = M[i] * rho[i] + Egamma[i];

  // m...m+n-1
  for (int j = 0; j < n; ++j) lhs[m + j] = ETrho[j] + N[j] * gamma[j];
}

void CG_for_CR_scaling(const std::vector<double>& b, std::vector<double>& x,
                       const std::vector<double>& M,
                       const std::vector<double>& N,
                       const std::vector<int>& ptr,
                       const std::vector<int>& rows) {
  int m = M.size();
  int n = N.size();

  // initial residual
  std::vector<double> r = b;

  // initial approximation
  x.assign(m + n, 0.0);

  // direction
  std::vector<double> p = r;

  int iter{};
  std::vector<double> Ap(m + n);

  while (iter < 100) {
    mult_by_CR_matrix(p, Ap, M, N, ptr, rows);

    double norm_r = dotProd(r, r);
    double alpha = norm_r / dotProd(p, Ap);

    // x = x + alpha * p
    vectorAdd(x, p, alpha);

    // r = r - alpha * Ap;
    vectorAdd(r, Ap, -alpha);

    // exit test
    if (norm2(r) / norm2(b) < 1e-6) break;

    double beta = dotProd(r, r) / norm_r;

    // p = r + beta * p
    vectorAdd(p, r, 1.0, beta);

    ++iter;
  }
}

void CurtisReidScaling(const std::vector<int>& ptr,
                       const std::vector<int>& rows,
                       const std::vector<double>& val, std::vector<int>& rowexp,
                       std::vector<int>& colexp) {
  // Takes as input the CSC matrix A.
  // Computes Curtis-Reid scaling exponents for the matrix, using powers of 2.

  int n = colexp.size();
  int m = rowexp.size();

  // rhs for CG
  std::vector<double> rhs(m + n, 0.0);

  // log2 abs b_i plus sum abs log2 A_i:
  double* sumlogrow = &rhs.data()[0];

  // log2 abs c_j plus sum abs log2 A_:j
  double* sumlogcol = &rhs.data()[m];

  // number of entries in each row and column
  std::vector<double> row_entries(m, 0.0);
  std::vector<double> col_entries(n, 0.0);

  // log A_ij
  for (int col = 0; col < n; ++col) {
    for (int el = ptr[col]; el < ptr[col + 1]; ++el) {
      int row = rows[el];
      if (val[el] != 0.0) {
        double temp = log2(std::abs(val[el]));
        sumlogrow[row] += temp;
        sumlogcol[col] += temp;
        row_entries[row] += 1.0;
        col_entries[col] += 1.0;
      }
    }
  }

  // solve linear system with CG
  std::vector<double> exponents(m + n);
  CG_for_CR_scaling(rhs, exponents, row_entries, col_entries, ptr, rows);

  // unpack exponents into various components
  for (int i = 0; i < m; ++i) rowexp[i] = -std::round(exponents[i]);
  for (int j = 0; j < n; ++j) colexp[j] = -std::round(exponents[m + j]);
}
