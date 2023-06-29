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

const double kSolutionDiffTolerance = 1e-6;
