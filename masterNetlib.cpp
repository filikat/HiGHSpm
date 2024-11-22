#include <cassert>
#include <cstring>  // For strchr
#include <iomanip>
#include <iostream>
#include <regex>
#include <sstream>

#include "Highs.h"
#include "Ipm.h"
#include "io/Filereader.h"
#include "parallel/HighsParallel.h"

enum ArgC { kMinArgC = 1, kOptionNlaArg = 1, kOptionFormat, kMaxArgC };

int main(int argc, char** argv) {
  // run all Netlib collection

  if (argc < kMinArgC || argc > kMaxArgC) {
    std::cerr << "======= How to use: ./test nla_option format_option "
                 "=======\n";
    std::cerr << "nla_option     : 0 aug sys, 1 norm eq\n";
    std::cerr << "format_option  : 0 full, 1 hybrid packed, 2 hybrid hybrid, 3 "
                 "packed packed\n";
    return 1;
  }

  std::ifstream netlib_names("../../Netlib/netlib_names.txt");
  std::string pb_name;
  std::string netlibPath = "../../Netlib/data/";

  Clock clock;
  int converged{};
  int total_problems{};

  std::stringstream ss{};
  ss << std::setw(20) << "Pb name";
  ss << std::setw(12) << "Status";
  ss << std::setw(6) << "iter";
  ss << std::setw(12) << "Time" << '\n';

  highs::parallel::initialize_scheduler();

  while (getline(netlib_names, pb_name)) {
    ++total_problems;
    Clock clock0;
    Clock clock1;
    // ===================================================================================
    // READ PROBLEM
    // ===================================================================================

    // Read LP using Highs MPS read
    Highs highs;
    //  highs.setOptionValue("output_flag", false);
    std::string model_file = pb_name;
    model_file.insert(0, netlibPath);
    std::string model = extractModelName(model_file);
    HighsStatus status = highs.readModel(model_file);
    assert(status == HighsStatus::kOk);
    double read_time = clock1.stop();
    const bool presolve = true;
    HighsLp lp;
    double presolve_time = -1;
    if (presolve) {
      clock1.start();
      status = highs.presolve();
      assert(status == HighsStatus::kOk);
      lp = highs.getPresolvedLp();
      presolve_time = clock1.stop();
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
    clock1.start();
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
      else if (lp.row_lower_[i] <= -kHighsInf &&
               lp.row_upper_[i] >= kHighsInf) {
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
      if (lp.row_lower_[i] != lp.row_upper_[i] &&
          lp.row_lower_[i] > -kHighsInf && lp.row_upper_[i] < kHighsInf) {
        // If both row_lower_ and row_upper_ are finite and different, add
        // slack:
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
    double setup_time = clock1.stop();

    // ===================================================================================
    // LOAD AND SOLVE THE PROBLEM
    // ===================================================================================
    clock1.start();

    // create instance of IPM
    Ipm ipm{};

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
    if (options.format < kOptionFormatMin ||
        options.format > kOptionFormatMax) {
      std::cerr << "Illegal value of " << options.format
                << " for option_format: must be in [" << kOptionFormatMin
                << ", " << kOptionFormatMax << "]\n";
      return 1;
    }

    // extract problem name without mps
    std::regex rgx("(.+)\\.mps");
    std::smatch match;
    std::regex_search(pb_name, match, rgx);
    pb_name = match[1];

    // load the problem
    ipm.load(n, m, obj.data(), rhs.data(), lower.data(), upper.data(),
             colptr.data(), rowind.data(), values.data(), constraints.data(),
             pb_name, options);
    double load_time = clock1.stop();

    // solve LP
    clock1.start();
    Output out = ipm.solve();
    double optimize_time = clock1.stop();
    if (out.status == "Optimal") {
      ++converged;
    }

    double run_time = clock0.stop();

    ss << std::setw(20) << pb_name << ' ';
    ss << std::setw(12) << out.status << ' ';
    ss << std::setw(6) << out.iterations << ' ';
    ss << std::setw(12) << std::fixed << std::setprecision(3) << optimize_time
       << '\n';
  }

  netlib_names.close();

  std::cout << ss.str();
  std::cout << "Converged: " << converged << " / " << total_problems << '\n';
  std::cout << "Time: " << clock.stop() << '\n';

  return 0;
}
