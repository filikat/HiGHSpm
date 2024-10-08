#ifndef HFACTOR_SOLVER_H
#define HFACTOR_SOLVER_H

#include "Highs.h"
#include "LinearSolver.h"

class HFactorSolver : public LinearSolver {
 public:
  // HFactor data
  std::vector<HighsInt> basic_index_;
  HFactor factor_;

  // Functions
  int factorAS(const HighsSparseMatrix& highs_a,
               const std::vector<double>& scaling) override;
  int factorNE(const HighsSparseMatrix& highs_a,
               const std::vector<double>& scaling) override;
  int solveNE(const std::vector<double>& rhs,
              std::vector<double>& lhs) override;
  int solveAS(const std::vector<double>& rhs_x,
              const std::vector<double>& rhs_y, std::vector<double>& lhs_x,
              std::vector<double>& lhs_y) override;
  void clear() override;
};

#endif
