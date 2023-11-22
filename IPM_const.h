#ifndef IPM_CONST_H
#define IPM_CONST_H

enum OptionNla {
  kOptionNlaMin = -1,
  kOptionNlaChoose = kOptionNlaMin,
  kOptionNlaCg,           // 0
  kOptionNlaAugmented,    // 1
  kOptionNlaNewton,       // 2
  kOptionNlaAugmentedCg,  // 3
  kOptionNlaNewtonCg,     // 4
  kOptionNlaMax = kOptionNlaNewtonCg,
  kOptionNlaDefault = kOptionNlaChoose
};

enum OptionMetis {
  kOptionMetisMin = -1,
  kOptionMetisChoose = kOptionMetisMin,
  kOptionMetisOff,
  kOptionMetisOn,
  kOptionMetisMax = kOptionMetisOn,
  kOptionMetisDefault = kOptionMetisChoose
};

enum OptionMaxDenseCol {
  kOptionMaxDenseColMin = 0,
  kOptionMaxDenseColMax = 100,
  kOptionMaxDenseColDefault = 1
};

enum OptionPredCor {
  kOptionPredCorMin = 0,
  kOptionPredCorAvoid = kOptionPredCorMin,
  kOptionPredCorUse = 1,
  kOptionPredCorMax = kOptionPredCorUse,
  kOptionPredCorDefault = kOptionPredCorUse
};

const int kMaxIterations = 100;
const double kIpmTolerance = 1e-6;

const double kOptionDenseColToleranceMin = 0;
const double kOptionDenseColToleranceDefault = 0.25;
const double kOptionDenseColToleranceMax = 1.1;

const double kSolutionDiffTolerance = 1e-6;

const double kCgTolerance = 1e-12;
const int kCgIterationLimit = 10000;

const double kPrimalRegularization = 0;
const double kDualRegularization = 0;

// set default solver
// 1: ssids
// 2: ma86
// 3: qdldl
// 4: cholmod
// 5 : highs
const int kDefaultSolver = 2;

// how to compute schur complement:
// 0 - MA86 with single rhs
// 1 - MA86 with multiple rhs  <=== use this
// 2 - Hfactor
const int kSchurMethod = 1;

const double kMetisSchurRelativeThreshold = 0.02;
const int kMetisSchurAbsThreshold = 1000;

#endif