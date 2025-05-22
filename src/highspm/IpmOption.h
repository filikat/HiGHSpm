#ifndef HIGHSPM_IPM_OPTION_H
#define HIGHSPM_IPM_OPTION_H

#include "IpmConst.h"
#include "io/HighsIO.h"

namespace highspm {

enum OptionNla {
  kOptionNlaMin = 0,
  kOptionNlaAugmented = kOptionNlaMin,
  kOptionNlaNormEq,
  kOptionNlaChoose,
  kOptionNlaMax = kOptionNlaChoose,
  kOptionNlaDefault = kOptionNlaChoose
};

enum OptionCrossover {
  kOptionCrossoverMin = 0,
  kOptionCrossoverOff = kOptionCrossoverMin,
  kOptionCrossoverOn,
  kOptionCrossoverMax = kOptionCrossoverOn,
  kOptionCrossoverDefault = kOptionCrossoverOff
};

enum OptionParallel {
  kOptionParallelMin = 0,
  kOptionParallelOff = kOptionParallelMin,  // tree off     node off
  kOptionParallelOn,                        // tree on      node on
  kOptionParallelChoose,                    // tree choose  node choose
  kOptionParallelTreeOnly,                  // tree on      node off
  kOptionParallelNodeOnly,                  // tree off     node on
  kOptionParallelMax = kOptionParallelNodeOnly,
  kOptionParallelDefault = kOptionParallelChoose
};

struct Options {
  // Solver options
  OptionNla nla = kOptionNlaDefault;
  OptionCrossover crossover = kOptionCrossoverDefault;
  OptionParallel parallel = kOptionParallelDefault;

  // Ipm parameters
  Int max_iter = kMaxIterDefault;
  double feasibility_tol = kIpmTolDefault;
  double optimality_tol = kIpmTolDefault;
  double crossover_tol = kIpmTolDefault;
  bool refine_with_ipx = true;
  double time_limit = -1.0;

  // Logging
  const HighsLogOptions* log_options;
  bool display = true;
};

}  // namespace highspm

#endif