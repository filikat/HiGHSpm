#include "Direct.h"
#include "IPM_caller.h"
#include "Highs.h"
#include "lp_data/HighsLpUtils.h"
#include "util/HighsMatrixPic.h"
#include "util/HighsUtils.h"
#include <filesystem>

bool infNormDiffOk(const std::vector<double> x0, const std::vector<double> x1) {
  assert(x1.size() >= x0.size());
  double norm_diff = 0;
  for (HighsInt ix = 0; ix < HighsInt(x0.size()); ix++)
    norm_diff = std::max(std::fabs(x0[ix] - x1[ix]), norm_diff);
  return norm_diff < 1e-12;
}

void syntheticTheta(std::vector<double> &theta) {
  int dim = theta.size();
  const double theta_log10_lo = -4;
  const double theta_log10_hi = -theta_log10_lo;
  const double theta_log10_dl = theta_log10_hi - theta_log10_lo;
  HighsRandom random;
  for (int i = 0; i < dim; i++) {
    double theta_log10 = theta_log10_lo + random.fraction() * theta_log10_dl;
    theta[i] = std::pow(10, theta_log10);
  }  
}

void callNewtonSolve(ExperimentData &experiment_data,
                    const HighsSparseMatrix &highs_a,
                    const std::vector<double> &theta,
                    const std::vector<double> &rhs,
                    const std::vector<double> &exact_sol,
                    const IpmOptions& ipm_options
                    ) {
                
  const int x_dim = highs_a.num_col_;
  const int y_dim = highs_a.num_row_;
  experiment_data.model_num_col = x_dim;
  experiment_data.model_num_row = y_dim;
  std::vector<double> lhs(y_dim);
  IpmInvert invert;
  invert.decomposer_source = ipm_options.decomposer_source;
  int newton_status = newtonInvert(highs_a, theta, invert,
				   ipm_options.max_dense_col,
                                   ipm_options.dense_col_tolerance,
				   experiment_data, true,
				   ipm_options.decomposer_source);
  if (!newton_status) {
    newton_status =
      newtonSolve(highs_a, theta, rhs, lhs, invert, experiment_data,
		  ipm_options.decomposer_source);
    experiment_data.nla_time.total += experiment_data.nla_time.solve;
    double solution_error = 0;
    for (int ix = 0; ix < y_dim; ix++)
      solution_error =
        std::max(std::fabs(exact_sol[ix] - lhs[ix]), solution_error);
    experiment_data.solution_error = solution_error;
    experiment_data.condition = newtonCondition(highs_a, theta, invert,
						ipm_options.decomposer_source);
  }
  invert.clear();
}

int main(int argc, char** argv){
  IPM_caller ipm{};
  if (!ipm.readOptionsOk(argc, argv)) return 1;
  ipm.reportOptions();

  IpmOptions& ipm_options = ipm.options_;

  int x_dim;
  int y_dim;
  HighsSparseMatrix matrix;

  // read A matrix from file
  const std::string model_file = ipm_options.model_file;
  std::filesystem::path p(model_file);
  std::string stem = p.stem().string();

  //  std::filesystem::path dir("../../result/");
  std::filesystem::path dir("result/");
  std::filesystem::path full_path = dir / stem;
  std::string new_path = full_path.string() + std::to_string(ipm_options.decomposer_source) + ".csv";
  std::cout << "\nExperimenting with file " << model_file << "\n";
  std::cout.flush();

  Highs highs;
  //highs.setOptionValue("output_flag", false);
  HighsStatus status = highs.readModel(model_file);
  const bool presolve = ipm.options_.presolve == kHighsOnString;
  HighsLp lp;
  double presolve_time = -1;
  double start_time;
  if (presolve) {
    start_time = getWallTime();
    status = highs.presolve();
    assert(status == HighsStatus::kOk);
    lp = highs.getPresolvedLp();
    presolve_time = getWallTime() - start_time;
  } else {
    lp = highs.getLp();
    presolve_time = 0;
  }
  // Scale the LP
  HighsOptions highs_options;

  highs_options.simplex_scale_strategy = kSimplexScaleStrategyMaxValue015;
  const bool force_scaling = true;
  scaleLp(highs_options, lp, force_scaling);


  highs_options.simplex_scale_strategy = kSimplexScaleStrategyMaxValue015;
  const bool scale_lp = true;
  if (scale_lp) {
    // Possibly force scaling (attempt) for well-scaled LPs
    const bool force_scaling = true;
    scaleLp(highs_options, lp, force_scaling);
    if (lp.is_scaled_) {
      analyseVectorValues(nullptr, "Column scale factors", lp.num_col_, lp.scale_.col);
      analyseVectorValues(nullptr, "Row scale factors", lp.num_row_, lp.scale_.row);
    }    
  }
  //  writeLpMatrixPicToFile(highs_options, "LpMatrix", lp);
  //  analyseMatrixSparsity(highs_options.log_options, "Matrix", lp.a_matrix_.num_col_, lp.a_matrix_.num_row_, lp.a_matrix_.start_, lp.a_matrix_.index_);
  if (status == HighsStatus::kOk) {
    matrix = highs.getLp().a_matrix_;
    y_dim = matrix.num_row_;
    int nnz = matrix.numNz();
    for (int ix = 0; ix < y_dim; ix++) {
      matrix.start_.push_back(++nnz);
      matrix.index_.push_back(ix);
      matrix.value_.push_back(1);
    }
    matrix.num_col_ += y_dim;
    x_dim = matrix.num_col_;
  } else {
    printf("HiGHS fails to read %s\n", model_file.c_str());
    // Use a test matrix if there's no ml.mps
    x_dim = 4;
    y_dim = 2;
    matrix.num_row_ = y_dim;
    matrix.num_col_ = x_dim;
    matrix.format_ = MatrixFormat::kRowwise;
    matrix.start_ = {0, 3, 6};
    matrix.index_ = {0, 1, 2, 0, 1, 3};
    matrix.value_ = {1, 1, 1, 1, -1, 1};
  }
  matrix.ensureColwise();
  HighsRandom random;
  const bool unit_solution = true;
  const double theta_base = 0.25;//1;//
  const double theta_random_mu = 0;//1e-3;   // 1e2;
  std::vector<double> theta;
  const bool synthetic_theta = false;//true;//
  if (synthetic_theta) {
    theta.resize(x_dim);
    syntheticTheta(theta);
  } else {
    for (int ix = 0; ix < x_dim; ix++) {
      const double theta_value = theta_base + theta_random_mu * random.fraction();
      theta.push_back(theta_value);
    }
  }

  // Test solution of
  //
  // [-\Theta^{-1} A^T ][x_star] = [rhs_x]
  // [      A       0  ][y_star]   [rhs_y]
  //
  // First directly, and then by solving
  //
  // A\Theta.A^T y_star = rhs_y + A\Theta.rhs_x
  //
  // before substituting x_star = \Theta(A^Ty_star - rhs_x)

  std::vector<double> x_star(x_dim);
  for (int ix = 0; ix < x_dim; ix++)
    x_star[ix] = unit_solution ? 1 : random.fraction();

  std::vector<double> y_star(y_dim);
  for (int ix = 0; ix < y_dim; ix++)
    y_star[ix] = unit_solution ? 1 : random.fraction();

  // Form rhs_x = -\Theta^{-1}.x_star + A^T.y_star
  std::vector<double> at_y_star;
  matrix.productTranspose(at_y_star, y_star);
  std::vector<double> rhs_x = x_star;
  for (int ix = 0; ix < x_dim; ix++) {
    rhs_x[ix] /= -theta[ix];
    rhs_x[ix] += at_y_star[ix];
  }
  // Form rhs_y = A.x_star
  std::vector<double> rhs_y;
  matrix.product(rhs_y, x_star);

  const bool augmented_solve =
    ipm_options.nla == kOptionNlaAugmented ||
    ipm_options.nla == kOptionNlaAugmentedCg;
  const bool newton_solve = !augmented_solve ||
    ipm_options.nla == kOptionNlaCg;

  assert(augmented_solve != newton_solve);
  std::vector<ExperimentData> experiment_data_list;
  if (newton_solve) {
    // Now solve the Newton equation
    //
    // Form rhs_newton == rhs_y + A\Theta.rhs_x
    
    std::vector<double> theta_rhs_x = rhs_x;
    for (int ix = 0; ix < x_dim; ix++)
      theta_rhs_x[ix] *= theta[ix];
    std::vector<double> a_theta_rhs_x;
    matrix.product(a_theta_rhs_x, theta_rhs_x);
    std::vector<double> rhs_newton = rhs_y;
    for (int ix = 0; ix < y_dim; ix++)
      rhs_newton[ix] += a_theta_rhs_x[ix];

    ExperimentData experiment_data;
    experiment_data.reset();
    experiment_data.model_name = ipm_options.model_file;
    callNewtonSolve(experiment_data, matrix, theta, rhs_newton, y_star, ipm_options);
    std::cout << experiment_data << "\n";
    experiment_data_list.push_back(experiment_data);
  }
  
  if (augmented_solve && ipm_options.decomposer_source != kOptionDecomposerSourceCholmod) {
    // Solve the augmented system
    std::vector<double> lhs_x;
    std::vector<double> lhs_y;

    ExperimentData experiment_data;
    experiment_data.reset();
    experiment_data.model_name = ipm_options.model_file;

    experiment_data.model_num_col = x_dim;
    experiment_data.model_num_row = y_dim;
    IpmInvert invert;
    invert.decomposer_source = ipm_options.decomposer_source;

    int augmented_status =
        augmentedInvert(matrix, theta, invert, experiment_data, ipm_options.decomposer_source);
    if (!augmented_status) {
      augmentedSolve(matrix, theta, rhs_x, rhs_y, lhs_x, lhs_y, invert,
		     experiment_data, ipm_options.decomposer_source);
      experiment_data.nla_time.total += experiment_data.nla_time.solve;
      double solution_error = 0;
      for (int ix = 0; ix < x_dim; ix++)
	solution_error =
          std::max(std::fabs(x_star[ix] - lhs_x[ix]), solution_error);
      for (int ix = 0; ix < y_dim; ix++)
	solution_error =
          std::max(std::fabs(y_star[ix] - lhs_y[ix]), solution_error);
      experiment_data.solution_error = solution_error;
      experiment_data.condition = augmentedCondition(matrix, theta, invert, ipm_options.decomposer_source);
    }
    invert.clear();
    std::cout << experiment_data << "\n";
    experiment_data_list.push_back(experiment_data);
  }


  writeDataToCSV(experiment_data_list, new_path);
  return 0;
}
