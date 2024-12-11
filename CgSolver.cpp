#include "CgSolver.h"

void CgSolver::clear() {}

int CgSolver::factorAS(const HighsSparseMatrix& A,
                       const std::vector<double>& scaling) {
  printf("Cg does not solve augmented system\n");
  return kDecomposerStatusErrorFactorise;
}
int CgSolver::solveAS(const std::vector<double>& rhs_x,
                      const std::vector<double>& rhs_y,
                      std::vector<double>& lhs_x, std::vector<double>& lhs_y) {
  printf("Cg does not solve augmented system\n");
  return kDecomposerStatusErrorSolve;
}

int CgSolver::factorNE(const HighsSparseMatrix& A,
                       const std::vector<double>& scaling) {
  assert(!valid_);

  M.reset(A, scaling);
  P.reset(A, scaling);

  valid_ = true;
  return kDecomposerStatusOk;
}

int CgSolver::solveNE(const std::vector<double>& rhs,
                      std::vector<double>& lhs) {
  assert(valid_);
  int iter = Cg(&M, &P, rhs, lhs, 1e-4, 100);
  printf("Cg iter %d\n", iter);
  return kDecomposerStatusOk;
}
