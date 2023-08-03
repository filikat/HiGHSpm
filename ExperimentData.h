#ifndef EXPERIMENTDATA_H_
#define EXPERIMENTDATA_H_
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "Highs.h"

const int kDataNotSet = -1;
const int kSystemTypeAugmented = 1;
const int kSystemTypeNewton = 2;

class ExperimentData {
public:
  std::string decomposer;
  std::string model_name;
  int model_num_col;
  int model_num_row;
  int model_num_dense_col;
  double model_max_dense_col;
  double dense_col_tolerance;
  int use_num_dense_col;
  int system_type;
  int system_size;
  // bool is_A_positive_definite;
  double system_max_dense_col;
  int system_nnz;
  int nnz_L;
  double solution_error;
  std::pair<double, double> residual_error;
  double fill_in_factor;
  double condition;

  // time
  double time_taken;
  double form_time;
  double setup_time;
  double analysis_time;
  double factorization_time;
  double solve_time;

  // Theta
  double theta_min;
  double theta_geomean;
  double theta_max;
  int theta_num_small;
  int theta_num_medium;
  int theta_num_large;
  int theta_order0;
  std::vector<int> theta_order_k;

  void reset() {
    decomposer = "na";
    model_num_col = kDataNotSet;
    model_num_row = kDataNotSet;
    model_num_dense_col = kDataNotSet;
    model_max_dense_col = kDataNotSet;
    dense_col_tolerance = kDataNotSet;
    use_num_dense_col = kDataNotSet;
    system_type = kDataNotSet;
    system_size = kDataNotSet;
    system_max_dense_col = kDataNotSet;
    system_nnz = kDataNotSet;
    nnz_L = kDataNotSet;
    fill_in_factor = kDataNotSet;
    condition = kDataNotSet;
    solution_error = kDataNotSet;
    residual_error.first = kDataNotSet;
    residual_error.second = kDataNotSet;
    time_taken = kDataNotSet;
    form_time = kDataNotSet;
    setup_time = kDataNotSet;
    analysis_time = kDataNotSet;
    factorization_time = kDataNotSet;
    solve_time = kDataNotSet;
    theta_min = kDataNotSet;
    theta_geomean = kDataNotSet;
    theta_max = kDataNotSet;
    theta_num_small = kDataNotSet;
    theta_num_medium = kDataNotSet;
    theta_num_large = kDataNotSet;
    theta_order0 = kDataNotSet;
    theta_order_k.clear();
  }
  void fillIn_LL();
  void fillIn_LDL();
  void analyseTheta(const std::vector<double> &theta);
  void reportAnalyseTheta(const std::vector<double> &theta);
};

double getWallTime();

std::ostream &operator<<(std::ostream &os, const ExperimentData &data);
void writeDataToCSV(const std::vector<ExperimentData> &data,
                    const std::string &filename);
std::pair<double, double> residualErrorAugmented(
    const HighsSparseMatrix &A, const std::vector<double> &theta,
    const std::vector<double> &rhs_x, const std::vector<double> &rhs_y,
    std::vector<double> &lhs_x, std::vector<double> &lhs_y);
std::pair<double, double> residualErrorNewton(const HighsSparseMatrix &A,
                                              const std::vector<double> &theta,
                                              const std::vector<double> &rhs,
                                              const std::vector<double> &lhs);
#endif
