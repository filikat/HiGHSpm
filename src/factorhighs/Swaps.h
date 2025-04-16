#ifndef FACTORHIGHS_SWAPS_H
#define FACTORHIGHS_SWAPS_H

void permuteWithSwaps(double* x, const int* swaps, int n, bool reverse = false);

void swapCols(char uplo, int n, double* A, int lda, int i, int j, int* swaps,
              int* sign);

void applySwaps(const int* swaps, int nrow, int ncol, double* R);

#endif