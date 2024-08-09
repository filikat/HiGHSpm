#ifndef CHOLMOD_SOLVER_H
#define CHOLMOD_SOLVER_H

#include "LinearSolver.h"
#include "cholmod.h"
#include <algorithm>

class CholmodSolver : public LinearSolver {
public:
  // Cholmod data
  cholmod_common c;
  cholmod_triplet *T;
  cholmod_sparse *a;
  cholmod_factor *L;
  cholmod_dense *x;
  cholmod_dense *b;

  // Functions
  int FactorAS(const HighsSparseMatrix &highs_a,
               const std::vector<double> &theta) override;
  int FactorNE(const HighsSparseMatrix &highs_a,
               const std::vector<double> &theta) override;
  int SolveNE(const HighsSparseMatrix &highs_a,
              const std::vector<double> &theta, const std::vector<double> &rhs,
              std::vector<double> &lhs) override;
  int SolveAS(const HighsSparseMatrix &highs_a,
              const std::vector<double> &theta,
              const std::vector<double> &rhs_x,
              const std::vector<double> &rhs_y, std::vector<double> &lhs_x,
              std::vector<double> &lhs_y) override;
  void Clear() override;
};

#endif
