#include "NormalEquations.h"

NormalEquations::NormalEquations(const HighsSparseMatrix &input_highs_a,
                                 const std::vector<double> &input_scaling)
  highs_a{input_highs_a}, scaling{input_scaling} {}

void NormalEquations::Apply(const std::vector<double> &rhs,
                            std::vector<double> &lhs) const {

  //  std::vector<double> temp(A.cols(), 0.0);
  std::vector<double> temp(highs_a.num_col_, 0.0);

  // temp = A^T * rhs
  highs_a.product(1.0, rhs, temp, true);

  // temp = temp * theta
  VectorDivide(temp, scaling);

  // lhs = A * temp
  highs_a.product(1.0, temp, lhs);
}
