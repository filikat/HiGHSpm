#include <cassert>
#include <cstring>  // For strchr
#include <iomanip>
#include <iostream>
#include <regex>
#include <sstream>

#include "Highs.h"
#include "highspm/Ipm.h"
#include "io/Filereader.h"
#include "ipm/IpxWrapper.h"
#include "parallel/HighsParallel.h"

enum ArgC {
  kMinArgC = 2,
  kDirectoryPath = 1,
  kOptionNlaArg,
  kOptionFormat,
  kOptionCrossover,
  kMaxArgC
};

int main(int argc, char** argv) {
  // run collection of problems

  if (argc < kMinArgC || argc > kMaxArgC) {
    std::cerr << "======= How to use: ./test problems_path nla_option "
                 "format_option crossover_option =======\n";
    std::cerr << "problems_path    : path to directory, it must contain "
                 "names.txt and subdirectory data/\n";
    std::cerr << "nla_option       : 0 aug sys, 1 norm eq\n";
    std::cerr << "format_option    : 0 full, 1 hybrid packed, 2 hybrid hybrid, "
                 "3 packed packed\n";
    std::cerr << "crossover_option : 0 off, 1 on\n";
    return 1;
  }

  std::string directory_path = argv[kDirectoryPath];
  if (directory_path.back() != '/') directory_path += "/";

  std::ifstream problems_names(directory_path + "names.txt");
  std::string pb_name;
  std::string problems_path = directory_path + "data/";

  Clock clock;
  Int converged{};
  Int total_problems{};

  std::stringstream ss{};
  ss << std::setw(30) << "Pb name";
  ss << std::setw(12) << "Status";
  ss << std::setw(12) << "iter";
  ss << std::setw(12) << "Time" << '\n';

  highs::parallel::initialize_scheduler();

  while (getline(problems_names, pb_name)) {
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
    model_file.insert(0, problems_path);
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
    //  constraints[i] is : =, <, >
    // ===================================================================================

    clock1.start();
    Int n, m;
    std::vector<double> obj, rhs, lower, upper, Aval;
    std::vector<Int> Aptr, Aind;
    std::vector<char> constraints;
    double offset;

    fillInIpxData(lp, n, m, offset, obj, lower, upper, Aptr, Aind, Aval, rhs,
                  constraints);

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

    // extract problem name without mps
    std::regex rgx("(.+)\\.(mps|lp)");
    std::smatch match;
    std::regex_search(pb_name, match, rgx);
    pb_name = match[1];

    // load the problem
    ipm.load(n, m, obj.data(), rhs.data(), lower.data(), upper.data(),
             Aptr.data(), Aind.data(), Aval.data(), constraints.data(), offset,
             pb_name, options);
    double load_time = clock1.stop();

    // solve LP
    clock1.start();
    IpmStatus ipm_status = ipm.solve();
    double optimize_time = clock1.stop();

    if ((options.crossover && ipm_status == kIpmStatusBasic) ||
        (!options.crossover && ipm_status == kIpmStatusPDFeas))
      ++converged;

    double run_time = clock0.stop();

    std::string status_string;
    switch (ipm_status) {
      case kIpmStatusError:
        status_string = "Error";
        break;
      case kIpmStatusMaxIter:
        status_string = "Max iter";
        break;
      case kIpmStatusNoProgress:
        status_string = "No progress";
        break;
      case kIpmStatusPDFeas:
        status_string = "PD feas";
        break;
      case kIpmStatusBasic:
        status_string = "Basic";
        break;
      default:
        break;
    }

    ss << std::setw(30) << pb_name << ' ';
    ss << std::setw(12) << status_string << ' ';
    ss << std::setw(10) << ipm.getIter() << ' ';
    ss << std::setw(12) << std::fixed << std::setprecision(3) << optimize_time
       << '\n';
  }

  problems_names.close();

  std::cout << ss.str();
  std::cout << "Converged: " << converged << " / " << total_problems << '\n';
  std::cout << "Time: " << clock.stop() << '\n';

  return 0;
}