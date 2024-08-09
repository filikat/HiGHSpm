#ifndef MA97_SOLVER_H
#define MA97_SOLVER_H

#include "LinearSolver.h"
#include "hsl_wrapper.h"

class MA97Solver : public LinearSolver {
public:
  // MA97 data
  void *akeep;
  void *fkeep;
  ma97_control_d control;
  ma97_info_d info;

  // MC68 data
  mc68_control control_perm;
  mc68_info info_perm;
  std::vector<int> order;

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
