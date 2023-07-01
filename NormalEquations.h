#ifndef NORMALEQUATIONS_H
#define NORMALEQUATIONS_H

#include "SparseMatrix.h"
#include "VectorOperations.h"
#include <vector>

class NormalEquations {

  //  const SparseMatrix &A{};
  const HighsSparseMatrix &highs_a{};
  const std::vector<double> &scaling{};

public:
  // constructor
  NormalEquations(//const SparseMatrix &input_A,
		  const HighsSparseMatrix &highs_a,
                  const std::vector<double> &input_scaling);

  // apply matrix: lhs = A * Theta * A^T * rhs
  void Apply(const std::vector<double> &rhs, std::vector<double> &lhs) const;
};

#endif
