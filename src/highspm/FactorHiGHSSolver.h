#ifndef HIGHSPM_FACTORHIGHS_SOLVER_H
#define HIGHSPM_FACTORHIGHS_SOLVER_H

#include <algorithm>

#include "LinearSolver.h"
#include "auxiliary/IntConfig.h"
#include "factorhighs/FactorHiGHS.h"

namespace highspm {

class FactorHiGHSSolver : public LinearSolver {
  // symbolic factorization
  Symbolic S_;

  // numeric factorization
  Numeric N_;

  // keep track of whether as or ne is being factorized
  bool use_as_ = true;

  Int choose(const HighsSparseMatrix& A, Options& options);
  Int setNla(const HighsSparseMatrix& A, Options& options);
  void setParallel(const Options& options);

 public:
  FactorHiGHSSolver(const Options& options);

  // Override functions
  Int factorAS(const HighsSparseMatrix& A,
               const std::vector<double>& scaling) override;
  Int factorNE(const HighsSparseMatrix& A,
               const std::vector<double>& scaling) override;
  Int solveNE(const std::vector<double>& rhs,
              std::vector<double>& lhs) override;
  Int solveAS(const std::vector<double>& rhs_x,
              const std::vector<double>& rhs_y, std::vector<double>& lhs_x,
              std::vector<double>& lhs_y) override;
  Int setup(const HighsSparseMatrix& A, Options& options) override;
  void clear() override;
  void finalise() override;
  double flops() const override;
  double spops() const override;
  double nz() const override;
};

}  // namespace highspm

#endif
