#ifndef NUMERIC_H
#define NUMERIC_H

#include <memory>
#include <vector>

#include "SolveHandler.h"
#include "Symbolic.h"

class Numeric {
  // columns of factorization, stored by supernode
  std::vector<std::vector<double>> sn_columns_{};

  // swaps of columns for each supernode, ordered locally within a block
  std::vector<std::vector<int>> swaps_{};

  // information about 2x2 pivots
  std::vector<std::vector<double>> pivot_2x2_{};

  // symbolic object
  const Symbolic& S_;

  // object to handle solve phase in different formats
  std::unique_ptr<SolveHandler> SH_;

  // lower triangle of original matrix, permuted
  std::vector<int> rowsA_{};
  std::vector<int> ptrA_{};
  std::vector<double> valA_{};

  friend class Factorise;

 public:
  // dynamic regularization applied to the matrix
  std::vector<double> total_reg_{};

  Numeric(const Symbolic& S);

  // Full solve with refinement
  void solve(std::vector<double>& x) const;

  // Iterative refinement
  void refine(const std::vector<double>& rhs, std::vector<double>& x) const;
  std::vector<double> residual(const std::vector<double>& rhs,
                               const std::vector<double>& x) const;
  std::vector<double> residualQuad(const std::vector<double>& rhs,
                                   const std::vector<double>& x) const;
  double computeOmega(const std::vector<double>& b,
                      const std::vector<double>& x) const;

  void conditionNumber() const;
};

#endif
