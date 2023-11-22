#include "Lapack_wrapper.h"

int Lapack_wrapper_factor(std::vector<double>& a, std::vector<int>& ipiv) {
  // Call Lapack dsytrf to factor the dense symmetric matrix contained in a.
  // a contains (at least) the lower triangular part of the matrix (col-wise).
  // On return, a and ipiv contain data to be used by dsytrs to solve linear
  // systems.

  char uplo = 'L';
  int n = ipiv.size();
  double wkopt;
  int lwork = -1;
  int info;

  // first call to query size of work
  dsytrf(&uplo, &n, a.data(), &n, ipiv.data(), &wkopt, &lwork, &info);
  if (info) return info;

  lwork = (int)wkopt;
  std::vector<double> work(lwork);

  // second call to factorize
  dsytrf(&uplo, &n, a.data(), &n, ipiv.data(), work.data(), &lwork, &info);
  if (info) return info;

  return 0;
}

int Lapack_wrapper_solve(std::vector<double>& a, std::vector<int>& ipiv,
                         std::vector<double>& b) {
  // Call Lapack dsytrs to solve the dense linear system with rhs b.
  // a and ipiv contain data previously computed by dsytrf.

  char uplo = 'L';
  int n = ipiv.size();
  int nrhs = 1;
  int info;

  // solve linear system
  dsytrs(&uplo, &n, &nrhs, a.data(), &n, ipiv.data(), b.data(), &n, &info);
  if (info) return info;

  return 0;
}