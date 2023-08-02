#include "Highs.h"
#include "io/Filereader.h"
#include "IPM_caller.h"
#include <cassert>
#include <iostream>
#include <cstring>  // For strchr


enum ArgC {
  kMinArgC = 2,
  kModelFileArg = 1,
  kOptionNlaArg,
  kOptionPredCorArg,
  kOptionMaxDenseColArg,
  kOptionDenseColToleranceArg,
  kMaxArgC
};

int main(int argc, char **argv) {

  if (argc < kMinArgC || argc > kMaxArgC) {
    std::cerr << "======= How to use: ./ipm LP_name.mps(.gz) nla_option predcor_option max_dense_col_option dense_col_tolerance_option =======\n";
    return 1;
  }

  // ===================================================================================
  // READ PROBLEM
  // ===================================================================================

  // Read LP using Highs MPS read
  Highs highs;
  //  highs.setOptionValue("output_flag", false);
  std::string model_file = argv[kModelFileArg];
  std::string model = extractModelName(model_file);
  HighsStatus status = highs.readModel(model_file);
  assert(status == HighsStatus::kOk);
  const bool presolve = false;
  HighsLp lp;
  if (presolve) {
    status = highs.presolve();
    assert(status == HighsStatus::kOk);
    lp = highs.getPresolvedLp();
  } else {
    lp = highs.getLp();
  }

  // ===================================================================================
  // CHANGE FORMULATION
  // ===================================================================================
  // Input problem must be in the form
  //
  //  min   obj^T * x
  //  s.t.  A * x {<=,=,>=} rhs
  //        lower <= x <= upper
  //
  //  constraints[i] is : 0 for  =
  //                      1 for >=
  //                     -1 for <=
  // ===================================================================================

  // Make a local copy of LP data to be modified
  int n = lp.num_col_;
  int m = lp.num_row_;
  std::vector<double> obj(lp.col_cost_);
  std::vector<double> lower(lp.col_lower_);
  std::vector<double> upper(lp.col_upper_);
  std::vector<int> colptr(lp.a_matrix_.start_);
  std::vector<int> rowind(lp.a_matrix_.index_);
  std::vector<double> values(lp.a_matrix_.value_);

  // Prepare vectors for different formulation
  std::vector<double> rhs(lp.num_row_);
  std::vector<int> constraints(lp.num_row_);

  int num_slacks{};

  // Set up constraints and rhs based on row_lower_ and row_upper_ for each
  // constraint
  for (int i = 0; i < m; ++i) {

    // equality constraint
    if (lp.row_lower_[i] == lp.row_upper_[i]) {
      constraints[i] = kConstraintTypeEqual;
      rhs[i] = lp.row_lower_[i];
    }

    // constraint <=
    else if (lp.row_lower_[i] <= -kHighsInf && lp.row_upper_[i] < kHighsInf) {
      constraints[i] = kConstraintTypeUpper;
      rhs[i] = lp.row_upper_[i];
    }

    // constraint >=
    else if (lp.row_lower_[i] > -kHighsInf && lp.row_upper_[i] >= kHighsInf) {
      constraints[i] = kConstraintTypeLower;
      rhs[i] = lp.row_lower_[i];
    }

    // no constraint
    else if (lp.row_lower_[i] <= -kHighsInf && lp.row_upper_[i] >= kHighsInf) {
      std::cout << "======= Free variable not yet supported =======\n";
      return 1;
    }

    // general constraint
    else {
      // keep track of how many slacks are needed to allocate memory
      ++num_slacks;
    }
  }

  // Reserve memory for slacks
  obj.reserve(n + num_slacks);
  lower.reserve(n + num_slacks);
  upper.reserve(n + num_slacks);
  colptr.reserve(n + num_slacks + 1);
  rowind.reserve(lp.a_matrix_.numNz() + num_slacks);
  values.reserve(lp.a_matrix_.numNz() + num_slacks);

  for (int i = 0; i < m; ++i) {

    if (lp.row_lower_[i] != lp.row_upper_[i] && lp.row_lower_[i] > -kHighsInf &&
        lp.row_upper_[i] < kHighsInf) {
      // If both row_lower_ and row_upper_ are finite and different, add slack:
      //   lb <= a^T * x <= ub
      //   becomes
      //   a^T * x - s = 0 and lb <= s <= ub.
      //
      //  This requires:
      //   - updating obj, lower, upper
      //   - adding a column of -identity to A

      constraints[i] = kConstraintTypeEqual;
      rhs[i] = 0.0;

      // add slack
      ++n;
      obj.push_back(0.0);
      lower.push_back(lp.row_lower_[i]);
      upper.push_back(lp.row_upper_[i]);
      colptr.push_back(colptr.back() + 1);
      rowind.push_back(i);
      values.push_back(-1.0);
    }
  }

  // ===================================================================================
  // LOAD AND SOLVE THE PROBLEM
  // ===================================================================================

  // create instance of IPM
  IPM_caller ipm{};

  // ===================================================================================
  // Identify the option values and check their validity
  // ===================================================================================
  ipm.option_nla =
      argc > kOptionNlaArg ? atoi(argv[kOptionNlaArg]) : kOptionNlaDefault;
  if (ipm.option_nla < kOptionNlaMin || ipm.option_nla > kOptionNlaMax) {
    std::cerr << "Illegal value of " << ipm.option_nla
              << " for option_nla: must be in [" << kOptionNlaMin << ", "
              << kOptionNlaMax << "]\n";
    return 1;
  }

  ipm.option_predcor = 
    argc > kOptionPredCorArg ? atoi(argv[kOptionPredCorArg]) : kOptionPredCorDefault;
  if (ipm.option_predcor < kOptionPredCorMin || ipm.option_predcor > kOptionPredCorMax) {
    std::cerr << "Illegal value of " << ipm.option_predcor
              << " for option_predcor: must be in [" << kOptionPredCorMin << ", "
              << kOptionPredCorMax << "]\n";
    return 1;
  }

  ipm.option_max_dense_col =
      argc > kOptionMaxDenseColArg ? atoi(argv[kOptionMaxDenseColArg]) : kOptionMaxDenseColDefault;
  if (ipm.option_max_dense_col < kOptionMaxDenseColMin || ipm.option_max_dense_col > kOptionMaxDenseColMax) {
    std::cerr << "Illegal value of " << ipm.option_max_dense_col
              << " for option_max_dense_col: must be in [" << kOptionMaxDenseColMin << ", "
              << kOptionMaxDenseColMax << "]\n";
    return 1;
  }

  ipm.option_dense_col_tolerance =
      argc > kOptionDenseColToleranceArg ? atoi(argv[kOptionDenseColToleranceArg]) : kOptionDenseColToleranceDefault;
  if (ipm.option_dense_col_tolerance < kOptionDenseColToleranceMin ||
      ipm.option_dense_col_tolerance > kOptionDenseColToleranceMax) {
    std::cerr << "Illegal value of " << ipm.option_dense_col_tolerance
              << " for option_dense_col_tolerance: must be in [" << kOptionDenseColToleranceMin << ", "
              << kOptionDenseColToleranceMax << "]\n";
    return 1;
  }

  // load the problem
  ipm.Load(n, m, obj.data(), rhs.data(), lower.data(), upper.data(),
           colptr.data(), rowind.data(), values.data(), constraints.data());

  // solve LP
  ipm.Solve();

  if (!ipm.experiment_data_record.empty()) {
    ipm.experiment_data_record[0].model_name = model;
    ipm.experiment_data_record[0].model_num_col = lp.num_col_;
    ipm.experiment_data_record[0].model_num_row = lp.num_row_;
    std::string csv_file_name = model + "_direct.csv";
    writeDataToCSV(ipm.experiment_data_record, csv_file_name);
  }
  return 0;
}
