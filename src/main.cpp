#include <cassert>
#include <cstring>  // For strchr
#include <iostream>
#include <regex>

#include "Highs.h"
#include "highspm/Ipm.h"
#include "io/Filereader.h"
#include "ipm/IpxWrapper.h"
#include "parallel/HighsParallel.h"

enum ArgC {
  kMinArgC = 2,
  kModelFileArg = 1,
  kOptionNlaArg,
  kOptionFormatArg,
  kOptionCrossoverArg,
  kOptionParallelArg,
  kMaxArgC
};

int main(int argc, char** argv) {
  if (argc < kMinArgC || argc > kMaxArgC) {
    std::cerr << "======= How to use: ./ipm LP_name.mps(.gz) nla_option "
                 "format_option crossover_option parallel_option =======\n";
    std::cerr << "nla_option       : 0 aug sys, 1 norm eq, 2 choose\n";
    std::cerr << "format_option    : 0 full, 1 hybrid packed, 2 hybrid hybrid, "
                 "3 packed packed\n";
    std::cerr << "crossover_option : 0 off, 1 on\n";
    std::cerr << "parallel_option : 0 off, 1 on, 2 choose, 3 tree only, 4 node "
                 "only\n";
    return 1;
  }

  highspm::Clock clock0;
  highspm::Clock clock;

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
  highspm::Int n, m;
  std::vector<double> obj, rhs, lower, upper, Aval;
  std::vector<highspm::Int> Aptr, Aind;
  std::vector<char> constraints;
  double offset;

  fillInIpxData(lp, n, m, offset, obj, lower, upper, Aptr, Aind, Aval, rhs,
                constraints);

  double setup_time = clock.stop();

  // ===================================================================================
  // IDENTIFY OPTIONS
  // ===================================================================================

  highspm::Options options{};

  // option to choose normal equations or augmented system
  options.nla = argc > kOptionNlaArg
                    ? (highspm::OptionNla)atoi(argv[kOptionNlaArg])
                    : highspm::kOptionNlaDefault;
  if (options.nla < highspm::kOptionNlaMin ||
      options.nla > highspm::kOptionNlaMax) {
    std::cerr << "Illegal value of " << options.nla
              << " for option_nla: must be in [" << highspm::kOptionNlaMin
              << ", " << highspm::kOptionNlaMax << "]\n";
    return 1;
  }

  // option to choose storage format inside FactorHiGHS
  options.format = argc > kOptionFormatArg
                       ? (highspm::OptionFormat)atoi(argv[kOptionFormatArg])
                       : highspm::kOptionFormatDefault;
  if (options.format < highspm::kOptionFormatMin ||
      options.format > highspm::kOptionFormatMax) {
    std::cerr << "Illegal value of " << options.format
              << " for option_format: must be in [" << highspm::kOptionFormatMin
              << ", " << highspm::kOptionFormatMax << "]\n";
    return 1;
  }

  // option to choose crossover
  options.crossover =
      argc > kOptionCrossoverArg
          ? (highspm::OptionCrossover)atoi(argv[kOptionCrossoverArg])
          : highspm::kOptionCrossoverDefault;
  if (options.crossover < highspm::kOptionCrossoverMin ||
      options.crossover > highspm::kOptionCrossoverMax) {
    std::cerr << "Illegal value of " << options.crossover
              << " for option_crossover: must be in ["
              << highspm::kOptionCrossoverMin << ", "
              << highspm::kOptionCrossoverMax << "]\n";
    return 1;
  }

  // option to choose parallel
  options.parallel =
      argc > kOptionParallelArg
          ? (highspm::OptionParallel)atoi(argv[kOptionParallelArg])
          : highspm::kOptionParallelDefault;
  if (options.parallel < highspm::kOptionParallelMin ||
      options.parallel > highspm::kOptionParallelMax) {
    std::cerr << "Illegal value of " << options.parallel
              << " for option_parallel: must be in ["
              << highspm::kOptionParallelMin << ", "
              << highspm::kOptionParallelMax << "]\n";
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
  highspm::Ipm ipm{};

  // scheduler should be already initialised, but no harm in doing it again
  // HighsTaskExecutor::shutdown(true);
  highs::parallel::initialize_scheduler();

  // load the problem
  ipm.load(n, m, obj.data(), rhs.data(), lower.data(), upper.data(),
           Aptr.data(), Aind.data(), Aval.data(), constraints.data(), offset,
           pb_name);
  double load_time = clock.stop();

  ipm.setOptions(options);

  // solve LP
  clock.start();
  ipm.solve();
  double optimize_time = clock.stop();

  highspm::IpmInfo info = ipm.getInfo();
  highspm::IpmStatus ipm_status = info.ipm_status;

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

  printf("\n");
  printf("Ipm iterations    : %d\n", info.ipm_iter);
  printf("Analyse NE time   : %.2f\n", info.analyse_NE_time);
  printf("Analyse AS time   : %.2f\n", info.analyse_AS_time);
  printf("Matrix time       : %.2f\n", info.matrix_time);
  printf("Factorisation time: %.2f\n", info.factor_time);
  printf("Solve time        : %.2f\n", info.solve_time);
  printf("Factorisations    : %d\n", info.factor_number);
  printf("Solves            : %d\n", info.solve_number);
  printf("Correctors        : %d\n", info.correctors);

  return 0;
}
