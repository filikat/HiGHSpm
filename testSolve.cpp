#include "Direct.h"
#include "Highs.h"
#include <boost/program_options.hpp>
#include <filesystem>
namespace po = boost::program_options;

bool infNormDiffOk(const std::vector<double> x0, const std::vector<double> x1) {
  assert(x1.size() >= x0.size());
  double norm_diff = 0;
  for (HighsInt ix = 0; ix < HighsInt(x0.size()); ix++)
    norm_diff = std::max(std::fabs(x0[ix] - x1[ix]), norm_diff);
  return norm_diff < 1e-12;
}

int callNewtonSolve(ExperimentData &experiment_data,
                    //const std::string model_file,
                    const HighsSparseMatrix &highs_a,
                    const std::vector<double> &theta,
                    const std::vector<double> &rhs,
                    const std::vector<double> &exact_sol,
                    const int option_max_dense_col,
                    const double option_dense_col_tolerance,
                    const int solver_type
                    ) {
                
  //ExperimentData experiment_data;
  const int x_dim = highs_a.num_col_;
  const int y_dim = highs_a.num_row_;
  std::vector<double> lhs(y_dim);
  IpmInvert invert;
  int newton_status = newtonInvert(highs_a, theta, invert, option_max_dense_col,
                                   option_dense_col_tolerance, experiment_data, true, solver_type);
  if (newton_status) return newton_status;
  newton_status =
      newtonSolve(highs_a, theta, rhs, lhs, invert, experiment_data, solver_type);
  if (newton_status) return newton_status;
  experiment_data.time_taken += experiment_data.solve_time;

  experiment_data.model_num_col = x_dim;
  experiment_data.model_num_row = y_dim;
  if (newton_status) {
    std::cout << experiment_data << "\n";
    invert.clear(solver_type);
    return newton_status;
  }
  double solution_error = 0;
  for (int ix = 0; ix < y_dim; ix++)
    solution_error =
        std::max(std::fabs(exact_sol[ix] - lhs[ix]), solution_error);
  experiment_data.solution_error = solution_error;
  experiment_data.condition = newtonCondition(highs_a, theta, invert, solver_type);
  //std::cout << experiment_data << "\n";
  return 0;
}

int main(int argc, char** argv){
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help", "produce help message")
    ("solver,s", po::value<int>(), "set the solver type")
    ("density,d", po::value<double>(), "set the density threshold")
    ("model,m", po::value<std::string>(), "model name") 
    ("max_dense_col,mdc", po::value<int>(), "set the maximum number of dense columns")
  ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);    

  if (vm.count("help")) {
    std::cout << desc << "\n";
    return 1;
  }

  int solverType = 2;  // default value
  if (vm.count("solver")) {
    solverType = vm["solver"].as<int>();
    std::cout << "Solver type was set to " << solverType << ".\n";
  } else {
    std::cout << "Solver type was not set. Using default: " << solverType << ".\n";
  }

  double densityThreshold = 0.5;  // default value
  if (vm.count("density")) {
    densityThreshold = vm["density"].as<double>();
    std::cout << "Density threshold was set to " << densityThreshold << ".\n";
  } else {
    std::cout << "Density threshold was not set. Using default: " << densityThreshold << ".\n";
  }

  std::string model;
  if (vm.count("model")) {
    model = vm["model"].as<std::string>();
  } else {
    std::cout << "Model was not set. Exiting...\n";
    return 1;
  }

  int maxDenseCol = 5;  // default value
  if (vm.count("max_dense_col")) {
    maxDenseCol = vm["max_dense_col"].as<int>();
    std::cout << "Maximum number of dense columns was set to " << maxDenseCol << ".\n";
  } else {
    std::cout << "Maximum number of dense columns was not set. Using default: " << maxDenseCol << ".\n";
  }

  int x_dim;
  int y_dim;
  HighsSparseMatrix matrix;
  Highs highs;
  //highs.setOptionValue("output_flag", false);

  //HighsStatus status = highs.readModel("ml.mps");
  // read A matrix from file
  const std::string model_file = model;
  std::filesystem::path p(model);
  std::string stem = p.stem().string();

  std::filesystem::path dir("../../result/");
  std::filesystem::path full_path = dir / stem;
  std::string new_path = full_path.string() + std::to_string(solverType)+".csv";
  HighsStatus status = highs.readModel(model_file);

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
  const bool unit_solution = false;
  double theta_random_mu = 1e-3;   // 1e2;
  std::vector<double> theta;
  for (int ix = 0; ix < x_dim; ix++) {
    const double theta_value = 1.0 + theta_random_mu * random.fraction();
    //    const double theta_value = 1.0 + std::pow(10, -(ix+1));
    theta.push_back(theta_value);
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

  const bool augmented_solve = true;
  const bool newton_solve = true;
  assert(augmented_solve || newton_solve);
  std::vector<ExperimentData> experiment_data_list;
  if (newton_solve) {
    ExperimentData experiment_data;
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

    //std::string new_path = full_path.string() + "_reduced.csv";
    int newton_status =
        callNewtonSolve(experiment_data, matrix, theta, rhs_newton, y_star, 0, 1.1, solverType);
    experiment_data_list.push_back(experiment_data);
    if (newton_status)
      return 1;
    /*
    newton_status = callNewtonSolve(matrix, theta, rhs_newton, y_star, 4, 0.4);
    if (newton_status) return 1;
    newton_status = callNewtonSolve(matrix, theta, rhs_newton, y_star, 10, 0.4);
    if (newton_status) return 1;
    newton_status = callNewtonSolve(matrix, theta, rhs_newton, y_star, 20, 0.4);
    if (newton_status) return 1;s
    */
    std::cout << experiment_data << "\n";

  }
  
  if (augmented_solve && solverType != 4) {
    // Solve the augmented system
    std::vector<double> lhs_x;
    std::vector<double> lhs_y;

    ExperimentData experiment_data;
    IpmInvert invert;

    int augmented_status =
        augmentedInvert(matrix, theta, invert, experiment_data, solverType);
    if (augmented_status) {
      invert.clear(solverType);
      return 1;
    }
    augmentedSolve(matrix, theta, rhs_x, rhs_y, lhs_x, lhs_y, invert,
                   experiment_data, solverType);
    experiment_data.time_taken += experiment_data.solve_time;
    experiment_data.model_num_col = x_dim;
    experiment_data.model_num_row = y_dim;
    if (augmented_status) {
      std::cout << experiment_data << "\n";
      return augmented_status;
    }
    double solution_error = 0;
    for (int ix = 0; ix < x_dim; ix++)
      solution_error =
          std::max(std::fabs(x_star[ix] - lhs_x[ix]), solution_error);
    for (int ix = 0; ix < y_dim; ix++)
      solution_error =
          std::max(std::fabs(y_star[ix] - lhs_y[ix]), solution_error);
    experiment_data.solution_error = solution_error;
    experiment_data.condition = augmentedCondition(matrix, theta, invert, solverType);
    std::cout << experiment_data << "\n";
    invert.clear();
    experiment_data_list.push_back(experiment_data);
    //std::string new_path = full_path.string() + "_augmented.csv";
    //writeDataToCSV(experiment_data_list, new_path);
  }


  writeDataToCSV(experiment_data_list, new_path);
  return 0;
}
