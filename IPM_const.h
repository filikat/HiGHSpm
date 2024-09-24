#ifndef IPM_CONST_H
#define IPM_CONST_H

enum OptionNla {
  kOptionNlaMin = 0,
  kOptionNlaAugmented = kOptionNlaMin,
  kOptionNlaNormEq,
  kOptionNlaMax = kOptionNlaNormEq,
  kOptionNlaDefault = kOptionNlaNormEq
};

enum OptionFact {
  kOptionFactMin = 0,
  kOptionFactChol = kOptionNlaMin,
  kOptionFactLDLt,
  kOptionFactMax = kOptionFactLDLt,
  kOptionFactDefault = kOptionFactLDLt,
};

enum OptionFormat {
  kOptionFormatMin = 0,
  kOptionFormatFull = kOptionFormatMin,
  kOptionFormatHybridPacked,
  kOptionFormatHybridHybrid,
  kOptionFormatMax = kOptionFormatHybridHybrid,
  kOptionFormatDefault = kOptionFormatHybridPacked
};

enum OptionPredCor {
  kOptionPredcorMin = 0,
  kOptionPredcorOff = kOptionPredcorMin,
  kOptionPredcorOn = 1,
  kOptionPredcorMax = kOptionPredcorOn,
  kOptionPredcorDefault = kOptionPredcorOn
};

enum Parameters { kParamNla, kParamFact, kParamFormat, kParamSize };

const int kMaxIterations = 100;
const double kIpmTolerance = 1e-8;

const double kSigmaInitial = 0.5;
const double kSigmaMin = 0.05;
const double kSigmaMax = 0.95;
const double kInteriorScaling = 0.999;

const double kPrimalStaticRegularization = 1e-12;
const double kDualStaticRegularization = 1e-8;

#endif