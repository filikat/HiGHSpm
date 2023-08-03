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
std::ostream &operator<<(std::ostream &os, const ExperimentData &data) {
  const int text_width = 30;
  const int num_width = 12;
  const int short_num_width = 4;
  const int pct_width = 4;
  const int int_pct_width = 3;
  const int two_dp = 2;
  double float_dim = double(data.system_size);
  assert(data.system_type != kDataNotSet);

  const double system_density =
      data.system_size ? 1e2 * double(data.system_nnz) / (float_dim * float_dim)
                       : -1;
  const double l_density =
      data.system_size && data.nnz_L >= 0
          ? 1e2 * double(data.nnz_L) /
                (float_dim * double(data.system_size + 1) * 0.5)
          : -1;
  const double sum_time = data.form_time + data.setup_time +
                          data.analysis_time + data.factorization_time +
                          data.solve_time;
  const double pct_sum_time =
      data.time_taken > 0 ? 1e2 * sum_time / data.time_taken : -1;
  const double pct_form_time =
      data.time_taken > 0 ? 1e2 * data.form_time / data.time_taken : -1;
  const double pct_setup_time =
      data.time_taken > 0 ? 1e2 * data.setup_time / data.time_taken : -1;
  const double pct_analysis_time =
      data.time_taken > 0 ? 1e2 * data.analysis_time / data.time_taken : -1;
  const double pct_factorization_time =
      data.time_taken > 0 ? 1e2 * data.factorization_time / data.time_taken
                          : -1;
  const double pct_solve_time =
      data.time_taken > 0 ? 1e2 * data.solve_time / data.time_taken : -1;
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
  os << std::left << std::setw(text_width) << "L nnz: " << std::right
     << std::setw(num_width) << data.nnz_L << " (" << std::right << std::fixed
     << l_density << "%)\n";
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
     << std::setw(num_width) << data.form_time << " (" << std::right
     << std::setw(int_pct_width) << roundDouble2Int(pct_form_time) << "%)\n"
     << std::left << std::setw(text_width) << "setup time: " << std::right
     << std::setw(num_width) << data.setup_time << " (" << std::right
     << std::setw(int_pct_width) << roundDouble2Int(pct_setup_time) << "%)\n"
     << std::left << std::setw(text_width) << "analyse time: " << std::right
     << std::setw(num_width) << data.analysis_time << " (" << std::right
     << std::setw(int_pct_width) << roundDouble2Int(pct_analysis_time) << "%)\n"
     << std::left << std::setw(text_width)
     << "factorization time: " << std::right << std::setw(num_width)
     << data.factorization_time << " (" << std::right
     << std::setw(int_pct_width) << roundDouble2Int(pct_factorization_time)
     << "%)\n"
     << std::left << std::setw(text_width) << "solve time: " << std::right
     << std::setw(num_width) << data.solve_time << " (" << std::right
     << std::setw(int_pct_width) << roundDouble2Int(pct_solve_time) << "%)\n"
     << std::left << std::setw(text_width) << "sum time: " << std::right
     << std::setw(num_width) << sum_time << " (" << std::right
     << std::setw(int_pct_width) << roundDouble2Int(pct_sum_time) << "%)\n"
     << std::left << std::setw(text_width) << "time taken: " << std::right
     << std::setw(num_width) << data.time_taken << "\n";
  return os;
}

void writeDataToCSV(const std::vector<ExperimentData> &data,
                    const std::string &filename) {
  if (data.empty())
    return;

  std::ofstream outputFile;
  outputFile.open(filename);

  const int system_type = data[0].system_type;
  outputFile << "Model," << data[0].model_name << "\n";
  outputFile << "Num col," << data[0].model_num_col << "\n";
  outputFile << "Num row," << data[0].model_num_row << "\n";
  if (system_type == kSystemTypeNewton) {
    outputFile << "max dense col," << data[0].model_max_dense_col << "\n";
    outputFile << "num dense col," << data[0].model_num_dense_col << "\n";
    outputFile << "dense col tolerance," << data[0].dense_col_tolerance << "\n";
  }
  outputFile << "System size," << data[0].system_size << "\n";

  // Write header
  outputFile << "Decomposer,";
  if (system_type == kSystemTypeNewton) {
    outputFile << "Num dense col,System max dense col,AAT NNZ,(%),";
  } else {
    outputFile << "System NNZ,(%),";
  }
  outputFile << "NNZ L,(%),Fill factor,Condition,Solution Error,Abs residual "
                "error,Rel residual error,";
  outputFile << "Time Taken, Form time, Setup time, Analyse time, "
                "Factorization time, Solve time\n";

  // Write data
  for (const auto &experimentData : data) {
    outputFile << experimentData.decomposer << ",";
    if (experimentData.system_type != system_type)
      break;
    double float_dim = double(experimentData.system_size);
    if (system_type == kSystemTypeNewton) {
      outputFile << experimentData.use_num_dense_col << ",";
      outputFile << experimentData.system_max_dense_col << ",";
    }
    outputFile << experimentData.system_nnz << ",";
    const double system_density =
        float_dim
            ? 1e2 * double(experimentData.system_nnz) / (float_dim * float_dim)
            : -1;
    outputFile << system_density << ",";

    outputFile << experimentData.nnz_L << ",";
    const double l_density =
        float_dim && experimentData.nnz_L >= 0
            ? 1e2 * double(experimentData.nnz_L) /
                  (double(float_dim) * double(experimentData.system_size + 1) *
                   0.5)
            : -1;
    outputFile << l_density << ",";
    outputFile << experimentData.fill_in_factor << ",";
    outputFile << experimentData.condition << ",";
    outputFile << experimentData.solution_error << ",";
    outputFile << experimentData.residual_error.first << ",";
    outputFile << experimentData.residual_error.second << ",";
    outputFile << experimentData.time_taken << ",";
    outputFile << experimentData.form_time << ",";
    outputFile << experimentData.setup_time << ",";
    outputFile << experimentData.analysis_time << ",";
    outputFile << experimentData.factorization_time << ",";
    outputFile << experimentData.solve_time << "\n";
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

void ExperimentData::fillIn_LL() {
  this->fill_in_factor =
      double(2 * this->nnz_L - this->system_size) / double(this->system_nnz);
}

void ExperimentData::fillIn_LDL() {
  this->fill_in_factor =
      double(2 * this->nnz_L + this->system_size) / double(this->system_nnz);
}

void ExperimentData::analyseTheta(const std::vector<double> &theta) {
  const int dim = theta.size();
  if (dim<=0) return;
  double min_log10_theta = kHighsInf;
  double max_log10_theta = 0;
  double theta_min = 0;
  double theta_geomean = 0;
  double theta_max = 0;
  for (int ix = 0; ix < dim; ix++) {
    theta_geomean += std::log10(theta[ix]);
    min_log10_theta = std::min(theta[ix], min_log10_theta);
    max_log10_theta = std::max(theta[ix], max_log10_theta);
  }
  assert(min_log10_theta>0);
  assert(max_log10_theta>0);
  this->theta_min = min_log10_theta;
  this->theta_max = max_log10_theta;
  min_log10_theta = std::log10(min_log10_theta);
  max_log10_theta = std::log10(max_log10_theta);
  double v;
  this->theta_order0 = int(std::floor(min_log10_theta));
  theta_geomean /= double(dim);
  theta_geomean = std::pow(10, theta_geomean);
  this->theta_num_small = 0;
  this->theta_num_medium = 0;
  this->theta_num_large = 0;
  for (int ix = 0; ix < dim; ix++) {
    if (theta[ix]*10 < theta_geomean) {
      this->theta_num_small++;
    } else if (theta[ix] <= theta_geomean*10) {
      this->theta_num_medium++;
    } else {
      this->theta_num_large++;
    }
  }
  this->theta_geomean = theta_geomean;
  int num_k = int(max_log10_theta) -  this->theta_order0 + 1;
  
  theta_order_k.assign(num_k, 0);
  for (int ix = 0; ix < dim; ix++) {
    int k = std::floor(std::log10(theta[ix]));
    k -= this->theta_order0;
    assert(k >= 0 && k < num_k);
    theta_order_k[k]++;
  }
  
}
void ExperimentData::reportAnalyseTheta(const std::vector<double> &theta) {
  this->analyseTheta(theta);
  const int dim = theta.size();
  printf("\nAnalysis of theta of dimension %d with values in [%g, %g] and geomean of %g\n",
	 dim, this->theta_min, this->theta_max, this->theta_geomean);
  printf("Number of [small, medium, large] theta values is [%d, %d, %d]\n",
	 int(this->theta_num_small), int(this->theta_num_medium), int(this->theta_num_large));
  const int num_k = theta_order_k.size();
  printf("Order |");
  for (int k = 0; k < num_k; k++) 
    printf(" %3d |", this->theta_order0+k);
  printf("\n");
  printf("Order |");
  for (int k = 0; k < num_k; k++) 
    printf(" %3d |", theta_order_k[k]);
  printf("\n");
}
  

