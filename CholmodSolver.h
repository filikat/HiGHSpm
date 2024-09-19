#ifndef CHOLMOD_SOLVER_H
#define CHOLMOD_SOLVER_H

#include <algorithm>

#include "LinearSolver.h"
#include "cholmod.h"

class CholmodSolver : public LinearSolver {
 public:
  // Cholmod data
  cholmod_common c_;
  cholmod_triplet* T_;
  cholmod_sparse* a_;
  cholmod_factor* L_;
  cholmod_dense* x_;
  cholmod_dense* b_;

  // Functions
  int factorAS(const HighsSparseMatrix& highs_a,
               const std::vector<double>& theta) override;
  int factorNE(const HighsSparseMatrix& highs_a,
               const std::vector<double>& theta) override;
  int solveNE(const HighsSparseMatrix& highs_a,
              const std::vector<double>& theta, const std::vector<double>& rhs,
              std::vector<double>& lhs) override;
  int solveAS(const HighsSparseMatrix& highs_a,
              const std::vector<double>& theta,
              const std::vector<double>& rhs_x,
              const std::vector<double>& rhs_y, std::vector<double>& lhs_x,
              std::vector<double>& lhs_y) override;
  void setup(const HighsSparseMatrix& A, const std::vector<int>& parameters) override{};
  void clear() override;
};

#endif
