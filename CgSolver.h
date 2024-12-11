#ifndef CG_SOLVER_H
#define CG_SOLVER_H

#include <algorithm>

#include "../FactorHiGHS/KrylovMethods.h"
#include "LinearSolver.h"

class CgSolver : public LinearSolver {
  IpmMatrix M;
  NeDiagPrec P;

 public:
  // Override functions
  int factorAS(const HighsSparseMatrix& A,
               const std::vector<double>& scaling) override;
  int factorNE(const HighsSparseMatrix& A,
               const std::vector<double>& scaling) override;
  int solveNE(const std::vector<double>& rhs,
              std::vector<double>& lhs) override;
  int solveAS(const std::vector<double>& rhs_x,
              const std::vector<double>& rhs_y, std::vector<double>& lhs_x,
              std::vector<double>& lhs_y) override;
  void clear() override;
};

#endif
