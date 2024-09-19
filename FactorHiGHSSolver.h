#ifndef FACTORHIGHS_SOLVER_H
#define FACTORHIGHS_SOLVER_H

#include <algorithm>

#include "../FactorHiGHS/Analyse.h"
#include "../FactorHiGHS/Factorise.h"
#include "LinearSolver.h"

class FactorHiGHSSolver : public LinearSolver {
 public:
  Symbolic S_;
  Numeric N_;

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
  void setup(const HighsSparseMatrix& A, const std::vector<int>& parameters) override;
  void clear() override;
};

#endif
