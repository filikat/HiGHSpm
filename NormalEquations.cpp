#include "NormalEquations.h"

NormalEquations::NormalEquations(const SparseMatrix &input_A,
				 const HighsSparseMatrix &input_highs_a,
                                 const std::vector<double> &input_scaling)
    : A{input_A}, highs_a{input_highs_a}, scaling{input_scaling} {}

void NormalEquations::Apply(const std::vector<double> &rhs,
                            std::vector<double> &lhs) const {

  std::vector<double> temp(A.cols(), 0.0);

  // temp = A^T * rhs
  mat_vec(A, highs_a, rhs, temp, 1.0, 't');

  // temp = temp * theta
  VectorDivide(temp, scaling);

  // lhs = A * temp
  mat_vec(A, highs_a, temp, lhs, 1.0, 'n');
}
