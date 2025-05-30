#include <cassert>
#include <cstring>  // For strchr
#include <iomanip>
#include <iostream>
#include <regex>
#include <sstream>

#include "Highs.h"
#include "highspm/HpmSolver.h"
#include "io/Filereader.h"
#include "ipm/IpxWrapper.h"
#include "parallel/HighsParallel.h"

enum ArgC {
  kMinArgC = 2,
  kDirectoryPath = 1,
  kOptionNlaArg,
  kOptionCrossoverArg,
  kOptionParallelArg,
  kMaxArgC
};

int main(int argc, char** argv) {
  // run collection of problems

  if (argc < kMinArgC || argc > kMaxArgC) {
    std::cerr << "======= How to use: ./test problems_path nla_option "
                 "crossover_option parallel_option =======\n";
    std::cerr << "problems_path    : path to directory, it must contain "
                 "names.txt and subdirectory data/\n";
    std::cerr << "nla_option       : 0 aug sys, 1 norm eq\n";
    std::cerr << "crossover_option : 0 off, 1 on\n";
    std::cerr << "parallel_option : 0 off, 1 on, 2 choose, 3 tree only, 4 node "
                 "only\n";
    return 1;
  }

  std::string directory_path = argv[kDirectoryPath];
  if (directory_path.back() != '/') directory_path += "/";

  std::ifstream problems_names(directory_path + "names.txt");
  std::string pb_name;
  std::string problems_path = directory_path + "data/";

  highspm::Clock clock;
  int converged{};
  int total_problems{};

  std::stringstream ss{};
  ss << std::setw(30) << "Pb name";
  ss << std::setw(12) << "Status";
  ss << std::setw(12) << "iter";
  ss << std::setw(12) << "Time" << '\n';

  highs::parallel::initialize_scheduler();

  while (getline(problems_names, pb_name)) {
    ++total_problems;
    highspm::Clock clock0;
    // ===================================================================================
    // READ PROBLEM
    // ===================================================================================

    // Read LP using Highs MPS read
    Highs highs;
    //  highs.setOptionValue("output_flag", false);
    std::string model_file = pb_name;
    model_file.insert(0, problems_path);
    std::string model = extractModelName(model_file);
    HighsStatus status = highs.readModel(model_file);
    if (status != HighsStatus::kOk) {
      printf("Problem reading file %s\n", model_file.c_str());
      continue;
    }

    const bool presolve = true;
    HighsLp lp;
    if (presolve) {
      status = highs.presolve();
      if (status != HighsStatus::kOk) {
        printf("Problem with presolve of file %s\n", model_file.c_str());
        continue;
      }
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
    // LOAD AND SOLVE THE PROBLEM
    // ===================================================================================

    highspm::HpmSolver hpm{};

    // ===================================================================================
    // Identify the option values and check their validity
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

    // Get information from highs object. I added some of these functions to
    // highs on purpose to extract the data. They will not be needed when hpm
    // will be built together with highs.
    const HighsOptions& hOptions = highs.getOptions();
    HighsCallback& hCallback = highs.getCallback();
    HighsTimer& hTimer = highs.getTimer();

    // start the highs timer, not needed once solver is built together with
    // highs.
    hTimer.start();

    // options.max_iter = 5;
    // options.refine_with_ipx = false;
    // options.time_limit = 1000;
    hpm.set(options, hOptions.log_options, hCallback, hTimer);

    // load the problem
    hpm.load(n, m, obj.data(), rhs.data(), lower.data(), upper.data(),
             Aptr.data(), Aind.data(), Aval.data(), constraints.data(), offset);

    hpm.solve();

    highspm::HpmInfo info = hpm.getInfo();
    highspm::IpmStatus ipm_status = info.ipm_status;

    if ((options.crossover && ipm_status == highspm::kIpmStatusBasic) ||
        (!options.crossover && ipm_status == highspm::kIpmStatusPDFeas))
      ++converged;

    double run_time = clock0.stop();

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

    std::string status_string;
    switch (ipm_status) {
      case highspm::kIpmStatusError:
        status_string = "Error";
        break;
      case highspm::kIpmStatusMaxIter:
        status_string = "Max iter";
        break;
      case highspm::kIpmStatusNoProgress:
        status_string = "No progress";
        break;
      case highspm::kIpmStatusPDFeas:
        status_string = "PD feas";
        break;
      case highspm::kIpmStatusBasic:
        status_string = "Basic";
        break;
      case highspm::kIpmStatusPrimalInfeasible:
        status_string = "Primal inf";
        break;
      case highspm::kIpmStatusDualInfeasible:
        status_string = "Dual inf";
        break;
      case highspm::kIpmStatusNotRun:
        status_string = "Not run";
        break;
      default:
        status_string = "Unknown";
        break;
    }

    ss << std::setw(30) << pb_name << ' ';
    ss << std::setw(12) << status_string << ' ';
    ss << std::setw(10) << info.ipm_iter << ' ';
    ss << std::setw(12) << std::fixed << std::setprecision(3) << run_time
       << '\n';
  }

  problems_names.close();

  std::cout << ss.str();
  std::cout << "Converged: " << converged << " / " << total_problems << '\n';
  std::cout << "Time: " << clock.stop() << '\n';

  return 0;
}