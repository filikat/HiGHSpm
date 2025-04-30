#ifndef HIGHSPM_IPM_OPTION_H
#define HIGHSPM_IPM_OPTION_H

#include "IpmConst.h"

namespace highspm {

enum OptionNla {
  kOptionNlaMin = 0,
  kOptionNlaAugmented = kOptionNlaMin,
  kOptionNlaNormEq,
  kOptionNlaChoose,
  kOptionNlaMax = kOptionNlaChoose,
  kOptionNlaDefault = kOptionNlaChoose
};

enum OptionFormat {
  kOptionFormatMin = 0,
  kOptionFormatFull = kOptionFormatMin,
  kOptionFormatHybridPacked,
  kOptionFormatHybridHybrid,
  kOptionFormatPackedPacked,
  kOptionFormatMax = kOptionFormatPackedPacked,
  kOptionFormatDefault = kOptionFormatHybridHybrid
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
  kOptionParallelChoose,                    // tree choose  node on
  kOptionParallelTreeOnly,                  // tree on      node off
  kOptionParallelNodeOnly,                  // tree off     node on
  kOptionParallelMax = kOptionParallelNodeOnly,
  kOptionParallelDefault = kOptionParallelChoose
};

struct Options {
  OptionNla nla = kOptionNlaDefault;
  OptionFormat format = kOptionFormatDefault;
  OptionCrossover crossover = kOptionCrossoverDefault;
  OptionParallel parallel = kOptionParallelDefault;
  Int max_iter = kMaxIterDefault;
  double feasibility_tol = kIpmTolDefault;
  double optimality_tol = kIpmTolDefault;
  double crossover_tol = kIpmTolDefault;
  bool refine_with_ipx = true;
  double time_limit = -1.0;
};

}  // namespace highspm

#endif