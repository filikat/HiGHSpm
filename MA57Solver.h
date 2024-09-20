#ifndef MA57_SOLVER_H
#define MA57_SOLVER_H

#include "LinearSolver.h"
#include "hsl_wrapper.h"

class MA57Solver : public LinearSolver {
 public:
  // MA57 data
  void* factors_;
  ma57_control_d control_;
  ma57_ainfo_d ainfo_;
  ma57_finfo_d finfo_;
  ma57_sinfo_d sinfo_;
  std::vector<int> order_;

  std::vector<int> col_;
  std::vector<int> row_;
  std::vector<double> val_;

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
  void clear() override;
};

#endif
