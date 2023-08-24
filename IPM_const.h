enum OptionDecomposerSource {
  kOptionDecomposerSourceMin = 0,
  kOptionDecomposerSourceCg = kOptionDecomposerSourceMin,
  kOptionDecomposerSourceSsids,   // 1
  kOptionDecomposerSourceHighs,      // 2
  kOptionDecomposerSourceMax = kOptionDecomposerSourceHighs,
  kOptionDecomposerSourceDefault = kOptionDecomposerSourceSsids
};

enum OptionNla {
  kOptionNlaMin = 0,
  kOptionNlaCg = kOptionNlaMin,
  kOptionNlaAugmented,   // 1
  kOptionNlaNewton,      // 2
  kOptionNlaAugmentedCg, // 3
  kOptionNlaNewtonCg,    // 4
  kOptionNlaMax = kOptionNlaNewtonCg,
  kOptionNlaDefault = kOptionNlaCg
};

enum OptionMaxDenseCol {
  kOptionMaxDenseColMin = 0,
  kOptionMaxDenseColMax = 100,
  kOptionMaxDenseColDefault = 1
};

enum DecomposerSource {
  kDecomposerSourceMin = 1,
  kDecomposerSourceSsids = kDecomposerSourceMin,
  kDecomposerSourceMa86,
  kDecomposerSourceQdldl,
  kDecomposerSourceCholmod,
  kDecomposerSourceHighs,
  kDecomposerSourceMax = kDecomposerSourceHighs
};

enum DecomposerStatus {
  kDecomposerStatusMin = 0,
  kDecomposerStatusOk = kDecomposerStatusMin,
  kDecomposerStatusErrorOom,
  kDecomposerStatusErrorFactorize,
  kDecomposerStatusErrorSolve,
  kDecomposerStatusErrorClear,
  kDecomposerStatusMax = kDecomposerStatusErrorClear
};

const double kOptionDenseColToleranceMin = 0;
const double kOptionDenseColToleranceDefault = 0.25;
const double kOptionDenseColToleranceMax = 1.1;

const double kSolutionDiffTolerance = 1e-6;

const double kCgTolerance = 1e-12;
const int kCgIterationLimit = 10000;

