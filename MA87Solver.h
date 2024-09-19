#ifndef MA87_SOLVER_H
#define MA87_SOLVER_H

#include "LinearSolver.h"
#include "hsl_wrapper.h"

class MA87Solver : public LinearSolver {
 public:
  // MA87 data
  void* keep_;
  ma87_control_d control_;
  ma87_info_d info_;

  // MC68 data
  mc68_control control_perm_;
  mc68_info info_perm_;
  std::vector<int> order_;

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
