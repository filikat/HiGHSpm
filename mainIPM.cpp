#include <cassert>
#include <cstring>  // For strchr
#include <iostream>
#include <regex>

#include "../FactorHiGHS/Auxiliary.h"
#include "Highs.h"
#include "Ipm.h"
#include "io/Filereader.h"
#include "parallel/HighsParallel.h"

enum ArgC {
  kMinArgC = 2,
  kModelFileArg = 1,
  kOptionNlaArg,
  kOptionFormat,
  kMaxArgC
};

int main(int argc, char** argv) {
  if (argc < kMinArgC || argc > kMaxArgC) {
    std::cerr << "======= How to use: ./ipm LP_name.mps(.gz) nla_option "
                 "format_option =======\n";
    std::cerr << "nla_option     : 0 aug sys, 1 norm eq\n";
    std::cerr << "format_option  : 0 full, 1 hybrid packed, 2 hybrid hybrid, 3 "
                 "packed packed\n";
    return 1;
  }

  Clock clock0;
  Clock clock;

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
  double read_time = clock.stop();
  const bool presolve = true;
  HighsLp lp;
  double presolve_time = -1;
  if (presolve) {
    clock.start();
    status = highs.presolve();
    assert(status == HighsStatus::kOk);
    lp = highs.getPresolvedLp();
    presolve_time = clock.stop();
  } else {
    lp = highs.getLp();
    presolve_time = 0;
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
  clock.start();
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

  int num_free_col = 0;
  for (int i = 0; i < n; ++i) {
    if (lower[i] <= -kHighsInf && upper[i] >= kHighsInf) num_free_col++;
  }
  if (num_free_col) {
    const double bound_on_free = 2e3;
    printf("Model has %d free columns\n",
           num_free_col);  //: replacing bounds with [%g, %g]\n",
                           // num_free_col, n, -bound_on_free, bound_on_free);
    /*for (int i = 0; i < n; ++i) {
      if (lower[i] <= -kHighsInf && upper[i] >= kHighsInf) {
        lower[i] = -bound_on_free;
        upper[i] = bound_on_free;
      }
    }*/
  }

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
  double setup_time = clock.stop();

  // ===================================================================================
  // LOAD AND SOLVE THE PROBLEM
  // ===================================================================================
  clock.start();

  // create instance of IPM
  Ipm ipm{};

  // scheduler should be already initialized, but no harm in doing it again
  // HighsTaskExecutor::shutdown(true);
  highs::parallel::initialize_scheduler();

  // ===================================================================================
  // Identify the option values and check their validity
  // ===================================================================================

  Options options{};

  // option to choose normal equations or augmented system
  options.nla =
      argc > kOptionNlaArg ? atoi(argv[kOptionNlaArg]) : kOptionNlaDefault;
  if (options.nla < kOptionNlaMin || options.nla > kOptionNlaMax) {
    std::cerr << "Illegal value of " << options.nla
              << " for option_nla: must be in [" << kOptionNlaMin << ", "
              << kOptionNlaMax << "]\n";
    return 1;
  }

  // option to choose storage format inside FactorHiGHS
  options.format =
      argc > kOptionFormat ? atoi(argv[kOptionFormat]) : kOptionFormatDefault;
  if (options.format < kOptionFormatMin || options.format > kOptionFormatMax) {
    std::cerr << "Illegal value of " << options.format
              << " for option_format: must be in [" << kOptionFormatMin << ", "
              << kOptionFormatMax << "]\n";
    return 1;
  }

  // extract problem name witout mps from path
  std::string pb_name{};
  std::regex rgx("([^/]+)\\.(mps|lp)");
  std::smatch match;
  std::regex_search(model_file, match, rgx);
  pb_name = match[1];

  // load the problem
  ipm.load(n, m, obj.data(), rhs.data(), lower.data(), upper.data(),
           colptr.data(), rowind.data(), values.data(), constraints.data(),
           pb_name, options);
  double load_time = clock.stop();

  // solve LP
  clock.start();
  Output output = ipm.solve();
  double optimize_time = clock.stop();

  double run_time = clock0.stop();

  if (run_time > 1e-3) {
    double sum_time =
        read_time + presolve_time + setup_time + load_time + optimize_time;
    printf("\nTime profile\n");
    printf("Read      %5.2f (%5.1f%% sum)\n", read_time,
           100 * read_time / sum_time);
    printf("Presolve  %5.2f (%5.1f%% sum)\n", presolve_time,
           100 * presolve_time / sum_time);
    printf("Setup     %5.2f (%5.1f%% sum)\n", setup_time,
           100 * setup_time / sum_time);
    printf("Load      %5.2f (%5.1f%% sum)\n", load_time,
           100 * load_time / sum_time);
    printf("Optimize  %5.2f (%5.1f%% sum)\n", optimize_time,
           100 * optimize_time / sum_time);
    printf("Sum       %5.2f (%5.1f%% run)\n", sum_time,
           100 * sum_time / run_time);
    printf("Run       %5.2f\n", run_time);
  }

  return 0;
}
