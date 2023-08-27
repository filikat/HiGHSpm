#include "lp_data/HighsLp.h"
#include "IPM_caller.h"
#include <cassert>
#include <cstring> // For strchr
#include <iostream>
//#include <fstream>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

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
  double start_time0 = getWallTime();
  // create instance of IPM
  IPM_caller ipm{};
  HighsLp lp;
  // ===================================================================================
  // READ OPTIONS AND PROBLEM 
  // ===================================================================================
  if (!ipm.readOptionsModelOk(argc, argv, lp)) return 1;
  ipm.reportOptions();

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
  double start_time = getWallTime();
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
    if (lower[i] <= -kHighsInf &&
	upper[i] >= kHighsInf) num_free_col++;
  }
  if (num_free_col) {
    const double bound_on_free = 2e3;
    printf("Model has %d/%d free columns: replacing bounds with [%g, %g]\n",
	   num_free_col, n, -bound_on_free, bound_on_free);
    for (int i = 0; i < n; ++i) {
      if (lower[i] <= -kHighsInf &&
	  upper[i] >= kHighsInf) {
	lower[i] = -bound_on_free;	
	upper[i] = bound_on_free;	
      }
    }
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
  double setup_time = getWallTime() - start_time;

  // ===================================================================================
  // LOAD AND SOLVE THE PROBLEM
  // ===================================================================================
  start_time = getWallTime();

  // load the problem
ipm.Load(n, m, obj.data(), rhs.data(), lower.data(), upper.data(),
           colptr.data(), rowind.data(), values.data(), constraints.data());
  double load_time = getWallTime() - start_time;

  // solve LP
  start_time = getWallTime();
  ipm.Solve();
  double optimize_time = getWallTime() - start_time;

  double run_time = getWallTime() - start_time0;
  double sum_time = ipm.read_time + ipm.presolve_time + setup_time + load_time + optimize_time;

  std::string csv_file_name = "result/" + ipm.options_.model_name +
    "_" + decomposerSource(ipm.options_.decomposer_source) +
    "_" + systemSolved(ipm.options_.nla);
  if (ipm.options_.nla == kOptionNlaNewton)
    csv_file_name += ("_" + std::to_string(ipm.options_.max_dense_col));
  csv_file_name += ".csv";
  std::ofstream output_stream;

  output_stream.open(csv_file_name);
  if (!output_stream)
    throw std::runtime_error("Cannot open file " + csv_file_name);

  
  if (!ipm.experiment_data_record.empty()) {
    ipm.experiment_data_record[0].model_name = ipm.options_.model_name;
    ipm.experiment_data_record[0].model_num_col = lp.num_col_;
    ipm.experiment_data_record[0].model_num_row = lp.num_row_;
    writeModelDataToCsv(ipm.experiment_data_record[0], output_stream);
    NlaTime sum_nla_time = sumNlaTime(ipm.experiment_data_record);
    double sum_sum_nla_time =
      sum_nla_time.form + sum_nla_time.setup +
      sum_nla_time.analysis + sum_nla_time.factorization +
      sum_nla_time.solve;
    if (sum_nla_time.total > 1e-3) {
      assert(sum_sum_nla_time>0);
      if (sum_sum_nla_time>0) {
	printf("NLA time profile\n");
	printf("Form      %5.2f (%5.1f%% sum)\n", sum_nla_time.form, 100*sum_nla_time.form/sum_sum_nla_time);
	printf("Setup     %5.2f (%5.1f%% sum)\n", sum_nla_time.setup, 100*sum_nla_time.setup/sum_sum_nla_time);
	printf("Analyse   %5.2f (%5.1f%% sum)\n", sum_nla_time.analysis, 100*sum_nla_time.analysis/sum_sum_nla_time);
	printf("Factorize %5.2f (%5.1f%% sum)\n", sum_nla_time.factorization, 100*sum_nla_time.factorization/sum_sum_nla_time);
	printf("Solve     %5.2f (%5.1f%% sum)\n", sum_nla_time.solve, 100*sum_nla_time.solve/sum_nla_time.total);
	printf("Total     %5.2f (%5.1f%% optimize)\n", sum_nla_time.total, 100*sum_nla_time.total/optimize_time);
      }
    }    

    if (sum_nla_time.total > 1e-3) {
      assert(sum_sum_nla_time>0);
      if (sum_sum_nla_time>0) {
	output_stream << "\nNLA time,Time,%\n";
	output_stream << "Form      ," << sum_nla_time.form << "," << 100*sum_nla_time.form/sum_sum_nla_time << "\n";
	output_stream << "Setup     ," << sum_nla_time.setup << "," << 100*sum_nla_time.setup/sum_sum_nla_time << "\n";
	output_stream << "Analyse   ," << sum_nla_time.analysis << "," << 100*sum_nla_time.analysis/sum_sum_nla_time << "\n";
	output_stream << "Factorize ," << sum_nla_time.factorization << "," << 100*sum_nla_time.factorization/sum_sum_nla_time << "\n";
	output_stream << "Solve     ," << sum_nla_time.solve << "," << 100*sum_nla_time.solve/sum_nla_time.total << "\n";
	output_stream << "Total     ," << sum_nla_time.total << "," << 100*sum_nla_time.total/optimize_time << "\n";
      }
    }    
    if (run_time > 1e-3) {
      output_stream << "\nRun time,Time,%\n";
      output_stream << "Read      ," << ipm.read_time << "," << 100*ipm.read_time/sum_time << "\n";
      output_stream << "Presolve  ," << ipm.presolve_time << "," << 100*ipm.presolve_time/sum_time << "\n";
      output_stream << "Setup     ," << setup_time << "," << 100*setup_time/sum_time << "\n";
      output_stream << "Load      ," << load_time << "," << 100*load_time/sum_time << "\n";
      output_stream << "Optimize  ," << optimize_time << "," << 100*optimize_time/sum_time << "\n";
      output_stream << "Sum       ," << sum_time << "," << 100*sum_time/run_time << "\n";
      output_stream << "Run       ," << run_time << "\n\n";
    }
    writeNlaDataToCsv(ipm.experiment_data_record, output_stream);
  }
  output_stream.close();
  if (run_time > 1e-3) {
    printf("\nTime profile\n");
    printf("Read      %5.2f (%5.1f%% sum)\n", ipm.read_time, 100*ipm.read_time/sum_time);
    printf("Presolve  %5.2f (%5.1f%% sum)\n", ipm.presolve_time, 100*ipm.presolve_time/sum_time);
    printf("Setup     %5.2f (%5.1f%% sum)\n", setup_time, 100*setup_time/sum_time);
    printf("Load      %5.2f (%5.1f%% sum)\n", load_time, 100*load_time/sum_time);
    printf("Optimize  %5.2f (%5.1f%% sum)\n", optimize_time, 100*optimize_time/sum_time);
    printf("Sum       %5.2f (%5.1f%% run)\n", sum_time, 100*sum_time/run_time);
    printf("Run       %5.2f\n", run_time);
  }

  return 0;
}
