#include <cassert>
#include <cstring>  // For strchr
#include <iostream>
#include <regex>

#include "Highs.h"
#include "highspm/HpmSolver.h"
#include "io/Filereader.h"
#include "ipm/IpxWrapper.h"
#include "parallel/HighsParallel.h"

enum ArgC {
  kMinArgC = 2,
  kModelFileArg = 1,
  kOptionNlaArg,
  kOptionCrossoverArg,
  kOptionParallelArg,
  kMaxArgC
};

int main(int argc, char** argv) {
  if (argc < kMinArgC || argc > kMaxArgC) {
    std::cerr << "======= How to use: ./ipm LP_name.mps(.gz) nla_option "
                 "crossover_option parallel_option =======\n";
    std::cerr << "nla_option       : 0 aug sys, 1 norm eq, 2 choose\n";
    std::cerr << "crossover_option : 0 off, 1 on\n";
    std::cerr << "parallel_option : 0 off, 1 on, 2 choose, 3 tree only, 4 node "
                 "only\n";
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
  const bool presolve = true;
  HighsLp lp;
  if (presolve) {
    status = highs.presolve();
    assert(status == HighsStatus::kOk);
    lp = highs.getPresolvedLp();
  } else {
    lp = highs.getLp();
  }

  // set highs options
  HighsOptions highs_options;
  highs_options.log_dev_level = 0;
  highs.passOptions(highs_options);

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

  highspm::Int n, m;
  std::vector<double> obj, rhs, lower, upper, Aval;
  std::vector<highspm::Int> Aptr, Aind;
  std::vector<char> constraints;
  double offset;

  fillInIpxData(lp, n, m, offset, obj, lower, upper, Aptr, Aind, Aval, rhs,
                constraints);

  // ===================================================================================
  // IDENTIFY OPTIONS
  // ===================================================================================

  highspm::HpmOptions options{};

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

  // Get information from highs object. I added some of these functions to highs
  // on purpose to extract the data. They will not be needed when hpm will be
  // built together with highs.
  const HighsOptions& hOptions = highs.getOptions();
  HighsCallback& hCallback = highs.getCallback();
  HighsTimer& hTimer = highs.getTimer();

  // start the highs timer, not needed once solver is built together with highs.
  hTimer.start();

  // ===================================================================================
  // LOAD AND SOLVE THE PROBLEM
  // ===================================================================================

  // create solver
  highspm::HpmSolver hpm{};

  // scheduler should be already initialised, but no harm in doing it again
  // HighsTaskExecutor::shutdown(true);
  highs::parallel::initialize_scheduler();

  // options.max_iter = 5;
  // options.refine_with_ipx = false;
  // options.time_limit = 1000;
  // options.display = 0;
  // options.timeless_log = true;
  hpm.set(options, hOptions.log_options, hCallback, hTimer);

  // load the problem
  hpm.load(n, m, obj.data(), rhs.data(), lower.data(), upper.data(),
           Aptr.data(), Aind.data(), Aval.data(), constraints.data(), offset);

  hpm.solve();

  highspm::HpmInfo info = hpm.getInfo();
  highspm::IpmStatus ipm_status = info.ipm_status;

  printf("\n");
  printf("Ipm iterations    : %" HIGHSINT_FORMAT "\n", info.ipm_iter);
  printf("Analyse NE time   : %.2f\n", info.analyse_NE_time);
  printf("Analyse AS time   : %.2f\n", info.analyse_AS_time);
  printf("Matrix time       : %.2f\n", info.matrix_time);
  printf("Factorisation time: %.2f\n", info.factor_time);
  printf("Solve time        : %.2f\n", info.solve_time);
  printf("Factorisations    : %" HIGHSINT_FORMAT "\n", info.factor_number);
  printf("Solves            : %" HIGHSINT_FORMAT "\n", info.solve_number);
  printf("Correctors        : %" HIGHSINT_FORMAT "\n", info.correctors);

  return 0;
}
