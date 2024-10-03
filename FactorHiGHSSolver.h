#ifndef FACTORHIGHS_SOLVER_H
#define FACTORHIGHS_SOLVER_H

#include <algorithm>

#include "../FactorHiGHS/Analyse.h"
#include "../FactorHiGHS/Factorise.h"
#include "LinearSolver.h"

class FactorHiGHSSolver : public LinearSolver {
 public:
  // symbolic and numeric factorization objects
  Symbolic S_;
  Numeric N_;

  // keep track of whether as or ne is being factorized
  bool use_as_ = true;

  // extreme values of the factorisation
  double maxD_{};
  double minD_{};
  double maxoffD_{};
  double minoffD_{};
  double worst_res_{};

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
  int setup(const HighsSparseMatrix& A,
            const std::vector<int>& parameters) override;
  void clear() override;
  void refine(const HighsSparseMatrix& A, const std::vector<double>& scaling,
              const std::vector<double>& rhs_x,
              const std::vector<double>& rhs_y, std::vector<double>& lhs_x,
              std::vector<double>& lhs_y) override;
  void finalise() override;

  void solveForRefineNE(const HighsSparseMatrix& A,
                        const std::vector<double>& scaling,
                        std::vector<double>& res_x, std::vector<double>& res_y);

  void extractData(LS_data& data) override;
};

#endif
