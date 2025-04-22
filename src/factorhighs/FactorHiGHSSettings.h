#ifndef FACTOR_HIGHS_SETTINGS_H
#define FACTOR_HIGHS_SETTINGS_H

#include <cmath>

#include "auxiliary/IntConfig.h"

// ===========================================================================
// SWITCHES
// ===========================================================================

#define PIVOTING
#define DATA_COLLECTION
// #define PRINT_ITER_REF
// #define PRINT_REGULARIZATION
// #define PRINT_CORRECTORS

// choose level of timing:
// - TIMING_0: no timing
// - TIMING_1: basic timing
// - TIMING_2: advanced timing
// - TIMING_3: extreme timing (timing of each BLAS call, considerably slower)
#define TIMING_1

// ===========================================================================
// PARAMETERS
// ===========================================================================

namespace highspm {

// supernode amalgamation
const Int kStartThreshRelax = 256;
const double kUpperRatioRelax = 0.02;
const double kLowerRatioRelax = 0.01;
const Int kMaxIterRelax = 10;
const Int kSnSizeRelax = 16;

// dense factorization
const Int kBlockSize = 128;
const double kAlphaBK = 0.01;  //(sqrt(17.0) + 1.0) / 8.0;
const Int kBlockGrainSize = 1;
const Int kBlockParallelThreshold = 5;

// regularization
const double kPrimalStaticRegularization = 1e-12;
const double kDualStaticRegularization = 1e-10;
const double kDynamicDiagCoeff = 1e-24;

// refinement
const Int kMaxRefinementIter = 3;
const double kRefinementTolerance = 1e-12;

}  // namespace highspm

#endif