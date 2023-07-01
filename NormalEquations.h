#ifndef NORMALEQUATIONS_H
#define NORMALEQUATIONS_H

#include "util/HighsSparseMatrix.h"
#include "VectorOperations.h"
#include <vector>

class NormalEquations {

  const HighsSparseMatrix &highs_a{};
  const std::vector<double> &scaling{};

public:
  // constructor
  NormalEquations(const HighsSparseMatrix &highs_a,
                  const std::vector<double> &input_scaling);

  // apply matrix: lhs = A * Theta * A^T * rhs
  void Apply(const std::vector<double> &rhs, std::vector<double> &lhs) const;
};

#endif
