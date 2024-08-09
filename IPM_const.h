#ifndef IPM_CONST_H
#define IPM_CONST_H

enum OptionNla {
  kOptionNlaMin = 0,
  kOptionNlaAugmented = kOptionNlaMin,
  kOptionNlaNewton, // 1
  kOptionNlaMax = kOptionNlaNewton,
  kOptionNlaDefault = kOptionNlaNewton
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

const double kPrimalRegularization = 0;
const double kDualRegularization = 0;

#endif