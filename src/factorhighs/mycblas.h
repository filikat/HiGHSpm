#ifndef FACTORHIGHS_MY_CBLAS_H
#define FACTORHIGHS_MY_CBLAS_H

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

void cblas_daxpy(const int n, const double alpha, const double* x,
                 const int incx, double* y, const int incy);
void cblas_dcopy(const int n, const double* x, const int incx, double* y,
                 const int incy);
void cblas_dscal(const int n, const double alpha, double* x, const int incx);
void cblas_dswap(const int n, double* x, const int incx, double* y,
                 const int incy);

// level 2

void cblas_dgemv(const enum CBLAS_ORDER order,
                 const enum CBLAS_TRANSPOSE transa, const int M, const int n,
                 const double alpha, const double* A, const int lda,
                 const double* x, const int incx, const double beta, double* y,
                 const int incy);

void cblas_dtpsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO uplo,
                 const enum CBLAS_TRANSPOSE transa, const enum CBLAS_DIAG diag,
                 const int n, const double* ap, double* x, const int incx);

void cblas_dtrsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO uplo,
                 const enum CBLAS_TRANSPOSE transa, const enum CBLAS_DIAG diag,
                 const int n, const double* a, const int lda, double* x,
                 const int incx);

void cblas_dger(const enum CBLAS_ORDER order, const int m, const int n,
                const double alpha, const double* x, const int incx,
                const double* y, const int incy, double* A, const int lda);

// level 3

void cblas_dgemm(const enum CBLAS_ORDER order,
                 const enum CBLAS_TRANSPOSE transa,
                 const enum CBLAS_TRANSPOSE transb, const int m, const int n,
                 const int k, const double alpha, const double* A,
                 const int lda, const double* B, const int ldb,
                 const double beta, double* C, const int ldc);

void cblas_dsyrk(const enum CBLAS_ORDER order, const enum CBLAS_UPLO uplo,
                 const enum CBLAS_TRANSPOSE trans, const int n, const int k,
                 const double alpha, const double* a, const int lda,
                 const double beta, double* C, const int ldc);

void cblas_dtrsm(const enum CBLAS_ORDER order, const enum CBLAS_SIDE side,
                 const enum CBLAS_UPLO uplo, const enum CBLAS_TRANSPOSE transa,
                 const enum CBLAS_DIAG diag, const int m, const int n,
                 const double alpha, const double* a, const int lda, double* b,
                 const int ldb);

#ifdef __cplusplus
}
#endif

#endif