#ifndef FACTORHIGHS_TIMING_H
#define FACTORHIGHS_TIMING_H

#include "FactorHiGHSSettings.h"

// defines for timing
#ifdef DEBUG

#if (defined(TIMING_1) || defined(TIMING_2) || defined(TIMING_3))
#define COARSE_TIMING
#endif

#if (defined(TIMING_2) || defined(TIMING_3))
#define FINE_TIMING
#endif

#if (defined(TIMING_3))
#define BLAS_TIMING
#endif

#endif

namespace highspm {

enum TimeItems {
  kTimeAnalyse,                           // TIMING_1
  kTimeAnalyseMetis,                      // TIMING_2
  kTimeAnalyseTree,                       // TIMING_2
  kTimeAnalyseCount,                      // TIMING_2
  kTimeAnalysePattern,                    // TIMING_2
  kTimeAnalyseSn,                         // TIMING_2
  kTimeAnalyseReorder,                    // TIMING_2
  kTimeAnalyseRelInd,                     // TIMING_2
  kTimeFactorise,                         // TIMING_1
  kTimeFactorisePrepare,                  // TIMING_2
  kTimeFactoriseAssembleOriginal,         // TIMING_2
  kTimeFactoriseAssembleChildrenFrontal,  // TIMING_2
  kTimeFactoriseAssembleChildrenClique,   // TIMING_2
  kTimeFactoriseDenseFact,                // TIMING_2
  kTimeDenseFact_main,                    // TIMING_2
  kTimeDenseFact_schur,                   // TIMING_2
  kTimeDenseFact_kernel,                  // TIMING_2
  kTimeDenseFact_convert,                 // TIMING_2
  kTimeDenseFact_pivoting,                // TIMING_2
  kTimeFactoriseTerminate,                // TIMING_2
  kTimeSolve,                             // TIMING_1
  kTimeSolvePrepare,                      // TIMING_2
  kTimeSolveSolve,                        // TIMING_2
  kTimeSolveResidual,                     // TIMING_2
  kTimeSolveOmega,                        // TIMING_2
  kTimeBlasStart,                         //
  kTimeBlas_copy = kTimeBlasStart,        // TIMING_3
  kTimeBlas_axpy,                         // TIMING_3
  kTimeBlas_scal,                         // TIMING_3
  kTimeBlas_swap,                         // TIMING_3
  kTimeBlas_gemv,                         // TIMING_3
  kTimeBlas_trsv,                         // TIMING_3
  kTimeBlas_tpsv,                         // TIMING_3
  kTimeBlas_ger,                          // TIMING_3
  kTimeBlas_trsm,                         // TIMING_3
  kTimeBlas_syrk,                         // TIMING_3
  kTimeBlas_gemm,                         // TIMING_3
  kTimeBlasEnd = kTimeBlas_gemm,          //
  kTimeSize                               //
};

}

#endif