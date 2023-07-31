enum OptionNla {
  kOptionNlaMin = 0,
  kOptionNlaCg = kOptionNlaMin,
  kOptionNlaAugmented,
  kOptionNlaNewton,
  kOptionNlaAugmentedCg,
  kOptionNlaNewtonCg,
  kOptionNlaMax = kOptionNlaNewtonCg,
  kOptionNlaDefault = kOptionNlaCg
};

enum OptionMaxDenseCol {
  kOptionMaxDenseColMin = 0,
  kOptionMaxDenseColMax = 100,
  kOptionMaxDenseColDefault = 0//1
};

enum OptionPredCor {
  kOptionPredCorMin = 0,
  kOptionPredCorAvoid = kOptionPredCorMin,
  kOptionPredCorUse = 1,
  kOptionPredCorMax = kOptionPredCorUse,
  kOptionPredCorDefault = kOptionPredCorAvoid
};

const double kOptionDenseColToleranceMin = 0;
const double kOptionDenseColToleranceDefault = 1.1;//0.5;
const double kOptionDenseColToleranceMax = 1.1;

const double kSolutionDiffTolerance = 1e-6;

const double kCgTolerance = 1e-12;
const int kCgIterationLimit = 10000;


