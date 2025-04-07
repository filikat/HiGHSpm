#ifndef PACKED_SOLVE_HANDLER_H
#define PACKED_SOLVE_HANDLER_H

#include "SolveHandler.h"

class PackedSolveHandler : public SolveHandler {
  void forwardSolve(std::vector<double>& x) const override;
  void backwardSolve(std::vector<double>& x) const override;
  void diagSolve(std::vector<double>& x) const override;

 public:
  PackedSolveHandler(const Symbolic& S,
                     const std::vector<std::vector<double>>& sn_columns);
};

#endif