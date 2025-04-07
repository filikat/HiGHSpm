#include "CallAndTimeBlas.h"
#include "DataCollector.h"

#include <cblas.h>

#include "Auxiliary.h"
#include "DenseFact.h"
#include "Timing.h"

// macros to interface with CBlas
#define TRANS(x) (x) == 'N' ? CblasNoTrans : CblasTrans
#define UPLO(x) (x) == 'U' ? CblasUpper : CblasLower
#define DIAG(x) (x) == 'N' ? CblasNonUnit : CblasUnit
#define SIDE(x) (x) == 'L' ? CblasLeft : CblasRight

// level 1

void callAndTime_daxpy(int n, double da, const double* dx, int incx, double* dy,
                       int incy) {
#ifdef BLAS_TIMING
  Clock clock;
#endif
  cblas_daxpy(n, da, dx, incx, dy, incy);
#ifdef BLAS_TIMING
  DataCollector::get()->sumTime(kTimeBlas_axpy, clock.stop());
#endif
}

void callAndTime_dcopy(int n, const double* dx, int incx, double* dy,
                       int incy) {
#ifdef BLAS_TIMING
  Clock clock;
#endif
  cblas_dcopy(n, dx, incx, dy, incy);
#ifdef BLAS_TIMING
  DataCollector::get()->sumTime(kTimeBlas_copy, clock.stop());
#endif
}

void callAndTime_dscal(int n, const double da, double* dx, int incx) {
#ifdef BLAS_TIMING
  Clock clock;
#endif
  cblas_dscal(n, da, dx, incx);
#ifdef BLAS_TIMING
  DataCollector::get()->sumTime(kTimeBlas_scal, clock.stop());
#endif
}

void callAndTime_dswap(int n, double* dx, int incx, double* dy, int incy) {
#ifdef BLAS_TIMING
  Clock clock;
#endif
  cblas_dswap(n, dx, incx, dy, incy);
#ifdef BLAS_TIMING
  DataCollector::get()->sumTime(kTimeBlas_swap, clock.stop());
#endif
}

// level 2

void callAndTime_dgemv(char trans, int m, int n, double alpha, const double* A,
                       int lda, const double* x, int incx, double beta,
                       double* y, int incy) {
#ifdef BLAS_TIMING
  Clock clock;
#endif
  cblas_dgemv(CblasColMajor, TRANS(trans), m, n, alpha, A, lda, x, incx, beta,
              y, incy);
#ifdef BLAS_TIMING
  DataCollector::get()->sumTime(kTimeBlas_gemv, clock.stop());
#endif
}

void callAndTime_dtpsv(char uplo, char trans, char diag, int n,
                       const double* ap, double* x, int incx) {
#ifdef BLAS_TIMING
  Clock clock;
#endif
  cblas_dtpsv(CblasColMajor, UPLO(uplo), TRANS(trans), DIAG(diag), n, ap, x,
              incx);
#ifdef BLAS_TIMING
  DataCollector::get()->sumTime(kTimeBlas_tpsv, clock.stop());
#endif
}

void callAndTime_dtrsv(char uplo, char trans, char diag, int n, const double* A,
                       int lda, double* x, int incx) {
#ifdef BLAS_TIMING
  Clock clock;
#endif
  cblas_dtrsv(CblasColMajor, UPLO(uplo), TRANS(trans), DIAG(diag), n, A, lda, x,
              incx);
#ifdef BLAS_TIMING
  DataCollector::get()->sumTime(kTimeBlas_trsv, clock.stop());
#endif
}

void callAndTime_dger(int m, int n, double alpha, const double* x, int incx,
                      const double* y, int incy, double* A, int lda) {
#ifdef BLAS_TIMING
  Clock clock;
#endif
  cblas_dger(CblasColMajor, m, n, alpha, x, incx, y, incy, A, lda);
#ifdef BLAS_TIMING
  DataCollector::get()->sumTime(kTimeBlas_ger, clock.stop());
#endif
}

// level 3

void callAndTime_dgemm(char transa, char transb, int m, int n, int k,
                       double alpha, const double* A, int lda, const double* B,
                       int ldb, double beta, double* C, int ldc) {
#ifdef BLAS_TIMING
  Clock clock;
#endif
  cblas_dgemm(CblasColMajor, TRANS(transa), TRANS(transb), m, n, k, alpha, A,
              lda, B, ldb, beta, C, ldc);
#ifdef BLAS_TIMING
  DataCollector::get()->sumTime(kTimeBlas_gemm, clock.stop());
#endif
}

void callAndTime_dsyrk(char uplo, char trans, int n, int k, double alpha,
                       const double* A, int lda, double beta, double* C,
                       int ldc) {
#ifdef BLAS_TIMING
  Clock clock;
#endif
  cblas_dsyrk(CblasColMajor, UPLO(uplo), TRANS(trans), n, k, alpha, A, lda,
              beta, C, ldc);
#ifdef BLAS_TIMING
  DataCollector::get()->sumTime(kTimeBlas_syrk, clock.stop());
#endif
}

void callAndTime_dtrsm(char side, char uplo, char trans, char diag, int m,
                       int n, double alpha, const double* A, int lda, double* B,
                       int ldb) {
#ifdef BLAS_TIMING
  Clock clock;
#endif
  cblas_dtrsm(CblasColMajor, SIDE(side), UPLO(uplo), TRANS(trans), DIAG(diag),
              m, n, alpha, A, lda, B, ldb);
#ifdef BLAS_TIMING
  DataCollector::get()->sumTime(kTimeBlas_trsm, clock.stop());
#endif
}
