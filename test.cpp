#include "Direct.h"
#include "Highs.h"

int main() {

  const int dim = 2;
  HighsSparseMatrix matrix;
  
  matrix.num_row_ = dim;
  matrix.num_col_ = 4;
  matrix.format_ = MatrixFormat::kRowwise;
  matrix.start_ = {0, 3, 6};
  matrix.index_ = {0, 1, 2, 0, 1, 3};
  matrix.value_ = {1, 1, 1, 1, -1, 1};
  const std::vector<double> theta = {1, 1, 1, 1};
  const HighsSparseMatrix AThetaAT;
  //=  computeAThetaAT_inner_product(matrix, theta.data());
  const std::vector<double> x_star = {1, 1};
  std::vector<double> rhs;
  rhs.assign(dim, 0);
  AThetaAT.product(1, x_star, rhs);
  std::vector<double> lhs(dim);
  //  newtonSolve(matrix, theta, rhs, lhs, 0, 0);
  double solution_error = 0;
  for (int ix = 0; ix < dim; ix++)
    solution_error = std::max(std::fabs(x_star[ix] - lhs[ix]), solution_error);
  std::vector<double> residual;
  residual = rhs;
  AThetaAT.product(-1, lhs, residual);
  double residual_error = 0;
  for (int ix = 0; ix < dim; ix++)
    residual_error = std::max(std::fabs(x_star[ix] - lhs[ix]), residual_error);
  assert(solution_error < 1e-6);
  assert(residual_error < 1e-6);
  
}

