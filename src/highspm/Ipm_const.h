#ifndef HIGHSPM_IPM_CONST_H
#define HIGHSPM_IPM_CONST_H

#include "auxiliary/IntConfig.h"

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

enum kOptionCrossover {
  kOptionCrossoverMin = 0,
  kOptionCrossoverOff = kOptionCrossoverMin,
  kOptionCrossoverOn,
  kOptionCrossoverMax = kOptionCrossoverOn,
  kOptionCrossoverDefault = kOptionCrossoverOff
};

enum kOptionParallel {
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
  Int nla = kOptionNlaDefault;
  Int format = kOptionFormatDefault;
  Int crossover = kOptionCrossoverDefault;
  Int parallel = kOptionParallelDefault;
};

enum IpmStatus {
  kIpmStatusError,
  kIpmStatusMaxIter,
  kIpmStatusNoProgress,
  kIpmStatusOptimal,
  kIpmStatusPDFeas,
  kIpmStatusBasic
};

enum LinearSolverStatus {
  kLinearSolverStatusOk = 0,
  kLinearSolverStatusErrorOom,
  kLinearSolverStatusErrorAnalyse,
  kLinearSolverStatusErrorFactorise,
  kLinearSolverStatusErrorSolve
};

// parameters for termination
const Int kMaxIterations = 100;
const double kIpmTolerance = 1e-8;
const Int kMaxBadIter = 5;

// parameters for correctors
const double kGammaCorrector = 0.1;
const double kSigmaAffine = 0.01;
const Int kMaxCorrectors = 5;
const double kMccIncreaseAlpha = 0.1;
const double kMccIncreaseMin = 0.1;
const double kSmallProduct = 1e-3;
const double kLargeProduct = 1e3;

// parameters for choice of AS or NE
const double kSpopsWeight = 30.0;
const double kRatioOpsThresh = 10.0;
const double kRatioSnThresh = 1.5;

// parameters for choice of parallelism
const double kLargeFlopsThresh = 1e7;
const double kLargeSpeedupThresh = 1.5;
const double kLargeSnThresh = 20.0;
const double kSmallSnThresh = 5.0;

}  // namespace highspm

#endif