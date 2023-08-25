#include "ExperimentData.h"
#include "Direct.h"
#include <iomanip>
//#include <ios>
double getWallTime() {
  using namespace std::chrono;
  using wall_clock = std::chrono::high_resolution_clock;
  return duration_cast<duration<double>>(wall_clock::now().time_since_epoch())
      .count();
}

int roundDouble2Int(const double value) { return int(value + 0.5); }

double systemDensity(const ExperimentData &data) {
  if (data.system_size <= 0 || data.system_nnz < 0) return -1;
  return 1e2 * double(data.system_nnz) / (double(data.system_size) * double(data.system_size));
}

double decompositionDensity(const ExperimentData &data) {
  if (data.system_size <= 0 || data.nnz_decomposition < 0) return -1;
  double density = -1;
  double full_decomposition = -1;
  if (data.decomposer == "HiGHS") {
    full_decomposition =
      double(data.system_size) * double(data.system_size);
  } else {
    full_decomposition = 
      0.5 * double(data.system_size) * double(data.system_size + 1);
  }
  density = 1e2 * double(data.nnz_decomposition) / full_decomposition;
  return density;
}
			    
std::ostream &operator<<(std::ostream &os, const ExperimentData &data) {
  const int text_width = 30;
  const int num_width = 12;
  const int short_num_width = 4;
  const int pct_width = 4;
  const int int_pct_width = 3;
  const int two_dp = 2;
  double float_dim = double(data.system_size);
  assert(data.system_type != kDataNotSet);

  const double system_density = systemDensity(data);
  const double decomposition_density = decompositionDensity(data);
  const double sum_time = data.nla_time.form + data.nla_time.setup +
                          data.nla_time.analysis + data.nla_time.factorization +
                          data.nla_time.solve;
  const double pct_sum_time =
      data.nla_time.total > 0 ? 1e2 * sum_time / data.nla_time.total : -1;
  const double pct_form_time =
      data.nla_time.total > 0 ? 1e2 * data.nla_time.form / data.nla_time.total : -1;
  const double pct_setup_time =
      data.nla_time.total > 0 ? 1e2 * data.nla_time.setup / data.nla_time.total : -1;
  const double pct_analysis_time =
      data.nla_time.total > 0 ? 1e2 * data.nla_time.analysis / data.nla_time.total : -1;
  const double pct_factorization_time =
      data.nla_time.total > 0 ? 1e2 * data.nla_time.factorization / data.nla_time.total
                          : -1;
  const double pct_solve_time =
      data.nla_time.total > 0 ? 1e2 * data.nla_time.solve / data.nla_time.total : -1;
  os << std::left << std::setw(text_width) << "model name:" << std::right
     << std::setw(num_width) << data.model_name << "\n"
     << std::left << std::setw(text_width) << "model num col:" << std::right
     << std::setw(num_width) << data.model_num_col << "\n"
     << std::left << std::setw(text_width) << "model num_row:" << std::right
     << std::setw(num_width) << data.model_num_row << "\n";
  if (data.system_type == kSystemTypeNewton) {
    os << std::left << std::setw(text_width) << "Newton system: ";
  } else {
    os << std::left << std::setw(text_width) << "Augmented system: ";
  }
  os << std::right << std::setw(num_width) << data.decomposer << "\n";

  os << std::fixed << std::setprecision(two_dp);
  if (data.system_type == kSystemTypeNewton) {
    os << std::left << std::setw(text_width)
       << "model max_dense_col:" << std::right << std::setw(num_width)
       << data.model_max_dense_col << "\n"
       << std::left << std::setw(23) << "model num_dense_col (@"
       << std::setw(short_num_width) << data.dense_col_tolerance
       << "): " << std::right << std::setw(num_width)
       << data.model_num_dense_col << "\n"
       << std::left << std::setw(text_width)
       << "use   num_dense_col:" << std::right << std::setw(num_width)
       << data.use_num_dense_col << "\n";
  }
  os << std::left << std::setw(text_width) << "system size: " << std::right
     << std::setw(num_width) << data.system_size << "\n";

  if (data.invert_status < 0) {
    os << std::left << std::setw(text_width) << "CG only\n";
    return os;
  } else if (data.invert_status > 0) {
    os << std::left << std::setw(text_width) << "Failure: "
       << std::right << std::setw(num_width) << data.invert_status << "\n";
    return os;
  }

  if (data.system_type == kSystemTypeNewton) {
    os << std::left << std::setw(text_width)
       << "system max col density: " << std::right << std::setw(num_width)
       << data.system_max_dense_col << "\n"
       << std::left << std::setw(text_width) << "AAT nnz: ";
  } else {
    os << std::left << std::setw(text_width) << "system nnz: ";
  }

  os << std::scientific;
  os << std::right << std::setw(num_width) << data.system_nnz << " ("
     << std::right << std::fixed << system_density << "%)\n";

  os << std::left << std::setw(text_width) << "Decomposition nnz: " << std::right
     << std::setw(num_width) << data.nnz_decomposition << " (" << std::right << std::fixed
     << decomposition_density << "%)\n";
  os << std::fixed;
  os << std::left << std::setw(text_width) << "fill-in: " << std::right
     << std::setw(num_width) << data.fill_in_factor << "\n";
  os << std::scientific;
  os << std::left << std::setw(text_width) << "condition: " << std::right
     << std::setw(num_width) << data.condition << "\n";
  os << std::left << std::setw(text_width) << "solution error: " << std::right
     << std::setw(num_width) << data.solution_error << "\n"
     << std::left << std::setw(text_width)
     << "rel (abs) residual error: " << std::right << std::setw(num_width)
     << data.residual_error.second << " (" << std::right
     << data.residual_error.first << ")\n";
  os << std::fixed;
  os << std::left << std::setw(text_width) << "form time: " << std::right
     << std::setw(num_width) << data.nla_time.form << " (" << std::right
     << std::setw(int_pct_width) << roundDouble2Int(pct_form_time) << "%)\n"
     << std::left << std::setw(text_width) << "setup time: " << std::right
     << std::setw(num_width) << data.nla_time.setup << " (" << std::right
     << std::setw(int_pct_width) << roundDouble2Int(pct_setup_time) << "%)\n"
     << std::left << std::setw(text_width) << "analyse time: " << std::right
     << std::setw(num_width) << data.nla_time.analysis << " (" << std::right
     << std::setw(int_pct_width) << roundDouble2Int(pct_analysis_time) << "%)\n"
     << std::left << std::setw(text_width)
     << "factorization time: " << std::right << std::setw(num_width)
     << data.nla_time.factorization << " (" << std::right
     << std::setw(int_pct_width) << roundDouble2Int(pct_factorization_time)
     << "%)\n"
     << std::left << std::setw(text_width) << "solve time: " << std::right
     << std::setw(num_width) << data.nla_time.solve << " (" << std::right
     << std::setw(int_pct_width) << roundDouble2Int(pct_solve_time) << "%)\n"
     << std::left << std::setw(text_width) << "sum time: " << std::right
     << std::setw(num_width) << sum_time << " (" << std::right
     << std::setw(int_pct_width) << roundDouble2Int(pct_sum_time) << "%)\n"
     << std::left << std::setw(text_width) << "time taken: " << std::right
     << std::setw(num_width) << data.nla_time.total << "\n";
  return os;
}

NlaTime sumNlaTime(const std::vector<ExperimentData> &experiment_data) {
  NlaTime sum_nla_time;
  sum_nla_time.form = 0;
  sum_nla_time.setup = 0;
  sum_nla_time.analysis = 0;
  sum_nla_time.factorization = 0;
  sum_nla_time.solve = 0;
  sum_nla_time.total = 0;
  for (const auto &data : experiment_data) {
    sum_nla_time.form += data.nla_time.form;
    sum_nla_time.setup += data.nla_time.setup;
    sum_nla_time.analysis += data.nla_time.analysis;
    sum_nla_time.factorization += data.nla_time.factorization;
    sum_nla_time.solve += data.nla_time.solve;
    sum_nla_time.total += data.nla_time.total;
  }
  return sum_nla_time; 
}

void writeDataToCSV(const std::vector<ExperimentData> &data,
                    const std::string &filename) {
  if (data.empty())
    return;

  std::ofstream outputFile;

  outputFile.open(filename);
  if (!outputFile)
  throw std::runtime_error("Cannot open file " + filename);

  const int system_type = data[0].system_type;
  outputFile << "Model," << data[0].model_name << "\n";
  outputFile << "Num col," << data[0].model_num_col << "\n";
  outputFile << "Num row," << data[0].model_num_row << "\n";
  if (system_type == kSystemTypeNewton) {
    outputFile << "max dense col," << data[0].model_max_dense_col << "\n";
    outputFile << "num dense col," << data[0].model_num_dense_col << "\n";
    outputFile << "dense col tolerance," << data[0].dense_col_tolerance << "\n";
  }

  // Write header
  outputFile << "Grep,Record,Decomposer,Model,System size,Min Theta,Max Theta,";
  if (system_type == kSystemTypeNewton) {
    outputFile << "Num dense col,System max dense col,Surplus large Theta,AAT NNZ,(%),";
  } else {
    outputFile << "System NNZ,(%),";
  }
  outputFile << "Decomposition NNZ,(%),Fill factor,Condition,Solution Error,Abs residual "
                "error,Rel residual error,";
  outputFile << "Time Taken, Form time, Setup time, Analyse time, "
                "Factorization time, Solve time\n";

  // Write data
  int record=0;
  for (const auto &experimentData : data) {
    record++;
    outputFile << "Grep,";
    outputFile << record << ",";
    outputFile << experimentData.decomposer << ",";
    outputFile << data[0].model_name << ",";
    outputFile << data[0].system_size << ",";
    outputFile << experimentData.theta_min << ",";
    outputFile << experimentData.theta_max << ",";
    if (experimentData.system_type != system_type) break;
    double float_dim = double(experimentData.system_size);
    if (system_type == kSystemTypeNewton) {
      outputFile << experimentData.use_num_dense_col << ",";
      outputFile << experimentData.system_max_dense_col << ",";
      outputFile << experimentData.theta_num_large - experimentData.system_size << ",";
    }

    if (experimentData.invert_status < 0) {
      outputFile << "CG only\n";
      continue;
    } else if (experimentData.invert_status > 0) {
      outputFile << "Failure," << experimentData.invert_status << "\n";
      continue;
    }

    outputFile << experimentData.system_nnz << ",";
    outputFile << systemDensity(experimentData) << ",";
    outputFile << experimentData.nnz_decomposition << ",";
    outputFile << decompositionDensity(experimentData) << ",";
    outputFile << experimentData.fill_in_factor << ",";
    outputFile << experimentData.condition << ",";
    outputFile << experimentData.solution_error << ",";
    outputFile << experimentData.residual_error.first << ",";
    outputFile << experimentData.residual_error.second << ",";
    outputFile << experimentData.nla_time.total << ",";
    outputFile << experimentData.nla_time.form << ",";
    outputFile << experimentData.nla_time.setup << ",";
    outputFile << experimentData.nla_time.analysis << ",";
    outputFile << experimentData.nla_time.factorization << ",";
    outputFile << experimentData.nla_time.solve << "\n";
  }

  outputFile.close();
}

std::pair<double, double> residualErrorAugmented(
    const HighsSparseMatrix &A, const std::vector<double> &theta,
    const std::vector<double> &rhs_x, const std::vector<double> &rhs_y,
    std::vector<double> &lhs_x, std::vector<double> &lhs_y) {
  std::vector<double> ATy;
  A.productTranspose(ATy, lhs_y);
  std::vector<double> Ax;
  A.product(Ax, lhs_x);
  std::pair<double, double> residual_error;
  residual_error.first = 0;
  double rhs_norm = 1; // so max(1, ||rhs||_\inf) is computed
  for (int ix = 0; ix < rhs_x.size(); ix++) {
    const double theta_i = !theta.empty() ? theta[ix] : 1;
    double residual = -1 / theta_i * lhs_x[ix] + ATy[ix] - rhs_x[ix];
    residual_error.first = std::max(std::fabs(residual), residual_error.first);
    rhs_norm = std::max(std::fabs(rhs_x[ix]), rhs_norm);
  }
  for (int ix = 0; ix < rhs_y.size(); ix++) {
    double residual = Ax[ix] - rhs_y[ix];
    residual_error.first = std::max(std::fabs(residual), residual_error.first);
    rhs_norm = std::max(std::fabs(rhs_y[ix]), rhs_norm);
  }
  residual_error.second = residual_error.first / rhs_norm;

  return residual_error;
}

std::pair<double, double> residualErrorNewton(const HighsSparseMatrix &A,
                                              const std::vector<double> &theta,
                                              const std::vector<double> &rhs,
                                              const std::vector<double> &lhs) {
  std::vector<double> AThetaATx;
  productAThetaAT(A, theta, lhs, AThetaATx);
  std::pair<double, double> residual_error;
  residual_error.first = 0;
  double rhs_norm = 1; // so max(1, ||rhs||_\inf) is computed
  for (int ix = 0; ix < rhs.size(); ix++) {
    residual_error.first =
        std::max(std::fabs(AThetaATx[ix] - rhs[ix]), residual_error.first);
    rhs_norm = std::max(std::fabs(rhs[ix]), rhs_norm);
  }
  residual_error.second = residual_error.first / rhs_norm;
  return residual_error;
}

void NlaTime::reset() {
  this->form = kDataNotSet;
  this->setup = kDataNotSet;
  this->analysis = kDataNotSet;
  this->factorization = kDataNotSet;
  this->solve = kDataNotSet;
  this->total = kDataNotSet;
}

void ExperimentData::fillIn_LL() {
  this->fill_in_factor = this->system_nnz ?
    double(2 * this->nnz_decomposition - this->system_size) / double(this->system_nnz) : -1;
}

void ExperimentData::fillIn_LU() {
  this->fill_in_factor = this->system_nnz ?
      double(this->nnz_decomposition) / double(this->system_nnz) : -1;
}

void ExperimentData::fillIn_LDL() {
  this->fill_in_factor = this->system_nnz ?
      double(2 * this->nnz_decomposition + this->system_size) / double(this->system_nnz) : -1;
}

void ExperimentData::analyseTheta(const std::vector<double> &theta, const bool quiet) {
  const int dim = theta.size();
  if (dim<=0) return;
  double min_log10_theta = kHighsInf;
  double max_log10_theta = 0;
  int theta_num_zero = 0;
  double theta_min = 0;
  double theta_geomean = 0;
  double theta_max = 0;
  for (int ix = 0; ix < dim; ix++) {
    if (theta[ix]) {
      theta_geomean += std::log10(theta[ix]);
      min_log10_theta = std::min(theta[ix], min_log10_theta);
      max_log10_theta = std::max(theta[ix], max_log10_theta);
    } else {
      theta_num_zero++;
    }
  }
  assert(min_log10_theta>0);
  assert(max_log10_theta>0);
  this->theta_min = min_log10_theta;
  this->theta_max = max_log10_theta;
  const double large_theta = this->theta_max * ok_theta_relative_tolerance;
  min_log10_theta = std::log10(min_log10_theta);
  max_log10_theta = std::log10(max_log10_theta);
  double v;
  this->theta_order0 = int(std::floor(min_log10_theta));
  theta_geomean /= double(dim-theta_num_zero);
  theta_geomean = std::pow(10, theta_geomean);
  this->theta_num_small = 0;
  this->theta_num_medium = 0;
  this->theta_num_large = 0;
  for (int ix = 0; ix < dim; ix++) {
    if (theta[ix] <= 0) continue;
    if (theta[ix]*10 < theta_geomean) {
      this->theta_num_small++;
    } else if (theta[ix] < large_theta) {
      this->theta_num_medium++;
    } else {
      this->theta_num_large++;
    }
  }
  this->theta_geomean = theta_geomean;
  int num_k = int(max_log10_theta) -  this->theta_order0 + 1;
  
  theta_order_k.assign(num_k, 0);
  for (int ix = 0; ix < dim; ix++) {
    if (theta[ix] <= 0) continue;
    int k = std::floor(std::log10(theta[ix]));
    k -= this->theta_order0;
    assert(k >= 0 && k < num_k);
    theta_order_k[k]++;
  }
  if (quiet) return;
  printf("\nAnalysis of theta of dimension %d with values in [%g, %g] and geomean of %g\n",
	 dim, this->theta_min, this->theta_max, this->theta_geomean);
  printf("Number of [zero, small, medium, large] theta values is [%d, %d, %d, %d]\n",
	 int(theta_num_zero), int(this->theta_num_small), int(this->theta_num_medium), int(this->theta_num_large));
  printf("Order |");
  for (int k = 0; k < num_k; k++) 
    if (theta_order_k[k])
      printf("%5d |", this->theta_order0+k);
  printf("\n");
  printf("Order |");
  for (int k = 0; k < num_k; k++) 
    if (theta_order_k[k])
      printf("%5d |", theta_order_k[k]);
  printf("\n");
}
