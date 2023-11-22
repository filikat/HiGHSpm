#ifndef LAPACK_WRAPPER_H
#define LAPACK_WRAPPER_H

#include <vector>

extern "C" void dsytrf(char* uplo, int* n, double* a, int* lda, int* ipiv,
                       double* work, int* lwork, int* info);
extern "C" void dsytrs(char* uplo, int* n, int* nrhs, double* a, int* lda,
                       int* ipiv, double* b, int* ldb, int* info);

int Lapack_wrapper_factor(std::vector<double>& a, std::vector<int>& ipiv);

int Lapack_wrapper_solve(std::vector<double>& a, std::vector<int>& ipiv,
                          std::vector<double>& b);

#endif