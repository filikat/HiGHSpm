#ifndef FACTORHIGHS_MY_CBLAS_H
#define FACTORHIGHS_MY_CBLAS_H

#include "auxiliary/IntConfig.h"

// Provide definition for cblas functions
// Based on Netlib implementation

enum CBLAS_ORDER { CblasRowMajor = 101, CblasColMajor = 102 };
enum CBLAS_TRANSPOSE {
  CblasNoTrans = 111,
  CblasTrans = 112,
  CblasConjTrans = 113
};
enum CBLAS_UPLO { CblasUpper = 121, CblasLower = 122 };
enum CBLAS_DIAG { CblasNonUnit = 131, CblasUnit = 132 };
enum CBLAS_SIDE { CblasLeft = 141, CblasRight = 142 };

#ifdef __cplusplus
extern "C" {
#endif

// level 1

void cblas_daxpy(const Int n, const double alpha, const double* x,
                 const Int incx, double* y, const Int incy);
void cblas_dcopy(const Int n, const double* x, const Int incx, double* y,
                 const Int incy);
void cblas_dscal(const Int n, const double alpha, double* x, const Int incx);
void cblas_dswap(const Int n, double* x, const Int incx, double* y,
                 const Int incy);

// level 2

void cblas_dgemv(const enum CBLAS_ORDER order,
                 const enum CBLAS_TRANSPOSE transa, const Int M, const Int n,
                 const double alpha, const double* A, const Int lda,
                 const double* x, const Int incx, const double beta, double* y,
                 const Int incy);

void cblas_dtpsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO uplo,
                 const enum CBLAS_TRANSPOSE transa, const enum CBLAS_DIAG diag,
                 const Int n, const double* ap, double* x, const Int incx);

void cblas_dtrsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO uplo,
                 const enum CBLAS_TRANSPOSE transa, const enum CBLAS_DIAG diag,
                 const Int n, const double* a, const Int lda, double* x,
                 const Int incx);

void cblas_dger(const enum CBLAS_ORDER order, const Int m, const Int n,
                const double alpha, const double* x, const Int incx,
                const double* y, const Int incy, double* A, const Int lda);

// level 3

void cblas_dgemm(const enum CBLAS_ORDER order,
                 const enum CBLAS_TRANSPOSE transa,
                 const enum CBLAS_TRANSPOSE transb, const Int m, const Int n,
                 const Int k, const double alpha, const double* A,
                 const Int lda, const double* B, const Int ldb,
                 const double beta, double* C, const Int ldc);

void cblas_dsyrk(const enum CBLAS_ORDER order, const enum CBLAS_UPLO uplo,
                 const enum CBLAS_TRANSPOSE trans, const Int n, const Int k,
                 const double alpha, const double* a, const Int lda,
                 const double beta, double* C, const Int ldc);

void cblas_dtrsm(const enum CBLAS_ORDER order, const enum CBLAS_SIDE side,
                 const enum CBLAS_UPLO uplo, const enum CBLAS_TRANSPOSE transa,
                 const enum CBLAS_DIAG diag, const Int m, const Int n,
                 const double alpha, const double* a, const Int lda, double* b,
                 const Int ldb);

#ifdef __cplusplus
}
#endif

#endif