#include "DgemmParallel.h"

#include "FactorHiGHSSettings.h"
#include "CallAndTimeBlas.h"
#include "parallel/HighsParallel.h"

dgemmParallelizer::dgemmParallelizer(const double* P, const double* R,
                                     double* Q, int col, int jb)
    : P_{P}, R_{R}, Q_{Q}, col_{col}, jb_{jb} {}

void dgemmParallelizer::run(int start, int end, double beta) const {
  callAndTime_dgemm('T', 'N', col_, end - start, jb_, -1.0, P_, jb_,
                    &R_[start * jb_], jb_, beta, &Q_[start * col_], col_);
}

void dgemmParallel(const double* P, const double* R, double* Q, int col, int jb,
                   int row, int nb, double beta) {
#ifdef PARALLEL_NODE
  // if there is enough work to be done, parallelize
  if (col >= nb / 2 && jb >= nb / 2 && row >= kBlockParallelThreshold * nb) {
    dgemmParallelizer gemmP(P, R, Q, col, jb);
    dgemmParallelizer* pt = &gemmP;

    // I need to use an object to call gemm, otherwise the task is too large and
    // static_assert in the parallel deque fails.
    highs::parallel::for_each(
        0, row, [pt, beta](int start, int end) { pt->run(start, end, beta); },
        kBlockGrainSize * jb);
  } else {
    callAndTime_dgemm('T', 'N', col, row, jb, -1.0, P, jb, R, jb, beta, Q, col);
  }
#else
  callAndTime_dgemm('T', 'N', col, row, jb, -1.0, P, jb, R, jb, beta, Q, col);
#endif
}