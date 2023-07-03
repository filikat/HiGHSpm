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
  kOptionMaxDenseColDefault = 1
};

const double kOptionDenseColToleranceMin = 0;
const double kOptionDenseColToleranceDefault = 0.5;
const double kOptionDenseColToleranceMax = 1.0;

const double kSolutionDiffTolerance = 1e-6;

const double kCgTolerance = 1e-12;
const int kCgInterationLimit = 5000;


