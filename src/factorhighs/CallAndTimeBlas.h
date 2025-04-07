#ifndef CALL_AND_TIME_BLAS_H
#define CALL_AND_TIME_BLAS_H

// level 1
void callAndTime_daxpy(int n, double da, const double* dx, int incx, double* dy,
                       int incy);
void callAndTime_dcopy(int n, const double* dx, int incx, double* dy, int incy);
void callAndTime_dscal(int n, const double da, double* dx, int incx);
void callAndTime_dswap(int n, double* dx, int incx, double* dy, int incy);

// level 2
void callAndTime_dgemv(char trans, int m, int n, double alpha, const double* A,
                       int lda, const double* x, int incx, double beta,
                       double* y, int incy);
void callAndTime_dtpsv(char uplo, char trans, char diag, int n,
                       const double* ap, double* x, int incx);
void callAndTime_dtrsv(char uplo, char trans, char diag, int n, const double* A,
                       int lda, double* x, int incx);
void callAndTime_dger(int m, int n, double alpha, const double* x, int incx,
                      const double* y, int incy, double* A, int lda);

// level 3
void callAndTime_dgemm(char transa, char transb, int m, int n, int k,
                       double alpha, const double* A, int lda, const double* B,
                       int ldb, double beta, double* C, int ldc);
void callAndTime_dsyrk(char uplo, char trans, int n, int k, double alpha,
                       const double* a, int lda, double beta, double* c,
                       int ldc);
void callAndTime_dtrsm(char side, char uplo, char trans, char diag, int m,
                       int n, double alpha, const double* a, int lda, double* b,
                       int ldb);

#endif