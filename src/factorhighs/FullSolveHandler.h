#ifndef FULL_SOLVE_HANDLER_H
#define FULL_SOLVE_HANDLER_H

#include "SolveHandler.h"

class FullSolveHandler : public SolveHandler {
  void forwardSolve(std::vector<double>& x) const override;
  void backwardSolve(std::vector<double>& x) const override;
  void diagSolve(std::vector<double>& x) const override;

 public:
  FullSolveHandler(const Symbolic& S,
                   const std::vector<std::vector<double>>& sn_columns);
};

#endif