#ifndef IPM_CONST_H
#define IPM_CONST_H

enum OptionNla {
  kOptionNlaMin = 0,
  kOptionNlaAugmented = kOptionNlaMin,
  kOptionNlaNormEq,
  kOptionNlaMax = kOptionNlaNormEq,
  kOptionNlaDefault = kOptionNlaNormEq
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

struct Options {
  int nla = kOptionNlaDefault;
  int format = kOptionFormatDefault;
};

// parameters for termination
const int kMaxIterations = 100;
const double kIpmTolerance = 1e-6;
const int kMaxBadIter = 5;

// parameters for correctors
const double kGammaCorrector = 0.1;
const double kSigmaAffine = 0.01;
const int kMaxCorrectors = 5;
const double kMccIncreaseAlpha = 0.1;
const double kMccIncreaseMin = 0.1;
const double kSmallProduct = 1e-3;
const double kLargeProduct = 1e3;

// other parameters
const double kInteriorScaling = 0.999;

#endif