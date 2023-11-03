#ifndef IPM_CONST_H
#define IPM_CONST_H

enum OptionNla {
  kOptionNlaMin = 0,
  kOptionNlaCg = kOptionNlaMin,
  kOptionNlaAugmented,       // 1
  kOptionNlaNewton,          // 2
  kOptionNlaAugmentedCg,     // 3
  kOptionNlaNewtonCg,        // 4
  kOptionNlaMetisAugmented,  // 5
  kOptionNlaMetisNormalEq,   // 6
  kOptionNlaMax = kOptionNlaMetisNormalEq,
  kOptionNlaDefault = kOptionNlaCg
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
  kOptionPredCorDefault = kOptionPredCorAvoid
};

const int kMaxIterations = 100;
const double kIpmTolerance = 1e-6;

const double kOptionDenseColToleranceMin = 0;
const double kOptionDenseColToleranceDefault = 0.25;
const double kOptionDenseColToleranceMax = 1.1;

const double kSolutionDiffTolerance = 1e-6;

const double kCgTolerance = 1e-12;
const int kCgIterationLimit = 10000;

const double kPrimalRegularization = 1e-6;
const double kDualRegularization = 1e-6;

// const double kPrimalRegularization = 0;
// const double kDualRegularization = 0;

#endif