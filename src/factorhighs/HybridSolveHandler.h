#ifndef HYBRID_SOLVE_HANDLER_H
#define HYBRID_SOLVE_HANDLER_H

#include "SolveHandler.h"

class HybridSolveHandler : public SolveHandler {
  const std::vector<std::vector<int>>& swaps_;
  const std::vector<std::vector<double>>& pivot_2x2_;

  void forwardSolve(std::vector<double>& x) const override;
  void backwardSolve(std::vector<double>& x) const override;
  void diagSolve(std::vector<double>& x) const override;

 public:
  HybridSolveHandler(const Symbolic& S,
                     const std::vector<std::vector<double>>& sn_columns,
                     const std::vector<std::vector<int>>& swaps,
                     const std::vector<std::vector<double>>& pivot_2x2);
};

#endif