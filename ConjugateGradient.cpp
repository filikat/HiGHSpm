#include "ConjugateGradient.h"
#include <cmath>

// Solve
//  A * lhs = rhs
// using Conjugate Gradient method without preconditioner.
// A is the normal equations matrix of an IPM.
void CG_solve(
    // INPUT
    const NormalEquations &A, const std::vector<double> &rhs, double tol,
    int maxit,
    // OUTPUT
    std::vector<double> &lhs) {

  int n = rhs.size();

  // initialize zero solution
  std::fill(lhs.begin(), lhs.end(), 0.0);

  std::vector<double> r(rhs);
  std::vector<double> p(r);
  double rho_old = DotProd(r, r);
  double rho_new{};
  double norm_rhs = Norm2(rhs);

  int iter{};
  while (iter < maxit) {
    std::vector<double> Ap(n, 0.0);
    A.Apply(p, Ap);
    double alpha = rho_old / DotProd(p, Ap);
    VectorAdd(lhs, p, alpha);
    VectorAdd(r, Ap, -alpha);
    rho_new = DotProd(r, r);
    if (std::sqrt(rho_new) < tol * norm_rhs) {
      break;
    }
    VectorAdd(p, r, rho_old / rho_new);
    VectorScale(p, rho_new / rho_old);
    rho_old = rho_new;
    ++iter;
  }

  // std::cout<<"CG iter: "<<iter<<'\n';
}