#include <cassert>
#include <cstring>  // For strchr
#include <iostream>
#include <regex>

#include "../FactorHiGHS/Auxiliary.h"
#include "Highs.h"
#include "Ipm.h"
#include "io/Filereader.h"
#include "ipm/IpxWrapper.h"
#include "parallel/HighsParallel.h"

enum ArgC {
  kMinArgC = 2,
  kModelFileArg = 1,
  kOptionNlaArg,
  kOptionFormat,
  kOptionCrossover,
  kMaxArgC
};

int main(int argc, char** argv) {
  if (argc < kMinArgC || argc > kMaxArgC) {
    std::cerr << "======= How to use: ./ipm LP_name.mps(.gz) nla_option "
                 "format_option crossover_option =======\n";
    std::cerr << "nla_option       : 0 aug sys, 1 norm eq\n";
    std::cerr << "format_option    : 0 full, 1 hybrid packed, 2 hybrid hybrid, "
                 "3 packed packed\n";
    std::cerr << "crossover_option : 0 off, 1 on\n";
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
  //  constraints[i] is : =, <, >
  // ===================================================================================

  clock.start();
  int n, m;
  std::vector<double> obj, rhs, lower, upper, Aval;
  std::vector<int> Aptr, Aind;
  std::vector<char> constraints;
  double offset;

  fillInIpxData(lp, n, m, offset, obj, lower, upper, Aptr, Aind, Aval, rhs,
                constraints);

  double setup_time = clock.stop();

  // ===================================================================================
  // IDENTIFY OPTIONS
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

  // option to choose crossover
  options.crossover = argc > kOptionCrossover ? atoi(argv[kOptionCrossover])
                                              : kOptionCrossoverDefault;
  if (options.crossover < kOptionCrossoverMin ||
      options.crossover > kOptionCrossoverMax) {
    std::cerr << "Illegal value of " << options.crossover
              << " for option_crossover: must be in [" << kOptionCrossoverMin
              << ", " << kOptionCrossoverMax << "]\n";
    return 1;
  }

  // extract problem name witout mps from path
  std::string pb_name{};
  std::regex rgx("([^/]+)\\.(mps|lp)");
  std::smatch match;
  std::regex_search(model_file, match, rgx);
  pb_name = match[1];

  // ===================================================================================
  // LOAD AND SOLVE THE PROBLEM
  // ===================================================================================
  clock.start();

  // create instance of IPM
  Ipm ipm{};

  // scheduler should be already initialized, but no harm in doing it again
  // HighsTaskExecutor::shutdown(true);
  highs::parallel::initialize_scheduler();

  // load the problem
  ipm.load(n, m, obj.data(), rhs.data(), lower.data(), upper.data(),
           Aptr.data(), Aind.data(), Aval.data(), constraints.data(), offset,
           pb_name, options);
  double load_time = clock.stop();

  // solve LP
  clock.start();
  IpmStatus ipm_status = ipm.solve();
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

  printf("Ipm iterations: %d\n", ipm.getIter());

  return 0;
}
