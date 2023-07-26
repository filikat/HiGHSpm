#include "Direct.h"
#include "Highs.h"

int main() {

  const bool use_lp = true;
  int dim;
  HighsSparseMatrix matrix;
  if (use_lp) {
    Highs highs;
    highs.setOptionValue("output_flag", false);
    HighsStatus status = highs.readModel("ml.mps");
    assert(status == HighsStatus::kOk);
    matrix = highs.getLp().a_matrix_;
    dim = matrix.num_row_;
    int nnz = matrix.numNz();
    for (int ix = 0; ix < dim; ix++) {
      matrix.start_.push_back(++nnz);
      matrix.index_.push_back(ix);
      matrix.value_.push_back(1);
    }
    matrix.num_col_ += dim;
  } else {
    dim = 2;
    matrix.num_row_ = dim;
    matrix.num_col_ = 4;
    matrix.format_ = MatrixFormat::kRowwise;
    matrix.start_ = {0, 3, 6};
    matrix.index_ = {0, 1, 2, 0, 1, 3};
    matrix.value_ = {1, 1, 1, 1, -1, 1};
  }
  HighsRandom random;
  std::vector<double> theta(matrix.num_col_);
  for (int ix = 0; ix < matrix.num_col_; ix++) theta[ix] = 1.0;// + 1e-3 * random.fraction();
  const HighsSparseMatrix AThetaAT = computeAThetaAT(matrix, theta);
  std::vector<double> x_star(dim);
  for (int ix = 0; ix < dim; ix++)
    x_star[ix] = random.fraction();

  std::vector<double> rhs;
  productAThetaAT(matrix, theta, x_star, rhs);
  
  std::vector<double> lhs(dim);
  ExperimentData data;
  newtonSolve(matrix, theta, rhs, lhs, 0, 0.8, data);

  double solution_error = 0;
  for (int ix = 0; ix < dim; ix++)
    solution_error = std::max(std::fabs(x_star[ix] - lhs[ix]), solution_error);
  double residual_error = residualErrorAThetaAT(matrix, theta, rhs, lhs);

  data.model_size = dim;
  data.solution_error = solution_error;
  data.residual_error = residual_error;
  std::cout << data << "\n";
  assert(solution_error < 1e-6);
  assert(residual_error < 1e-6);
}

