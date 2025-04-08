#include "CallAndTimeBlas.h"
#include "DataCollector.h"

void permuteWithSwaps(double* x, const int* swaps, int n, bool reverse) {
    // Apply swaps to vector x of length n
  
    if (!reverse) {
      // apply the swaps in forward order
      for (int i = 0; i < n; ++i) {
        if (swaps[i] != i) std::swap(x[i], x[swaps[i]]);
      }
    } else {
      // apply the swaps in backward order
      for (int i = n - 1; i >= 0; --i) {
        if (swaps[i] != i) std::swap(x[i], x[swaps[i]]);
      }
    }
  }

void swapCols(char uplo, int n, double* A, int lda, int i, int j, int* swaps,
              int* sign) {
  // Exchange rows/cols i and j of symmetric matrix A

  // make sure that i < j
  if (i == j) return;
  if (i > j) std::swap(i, j);

  // swap diagonal elements
  std::swap(A[i + i * lda], A[j + j * lda]);

  // swap rest of rows/cols
  if (uplo == 'L') {
    callAndTime_dswap(i, &A[i], lda, &A[j], lda);
    callAndTime_dswap(n - j - 1, &A[j + 1 + i * lda], 1, &A[j + 1 + j * lda],
                      1);
    callAndTime_dswap(j - i - 1, &A[i + 1 + i * lda], 1, &A[j + (i + 1) * lda],
                      lda);
  } else {
    callAndTime_dswap(i, &A[i * lda], 1, &A[j * lda], 1);
    callAndTime_dswap(n - j - 1, &A[i + (j + 1) * lda], lda,
                      &A[j + (j + 1) * lda], lda);
    callAndTime_dswap(j - i - 1, &A[i + (i + 1) * lda], lda,
                      &A[i + 1 + j * lda], 1);
  }

  // swap pivot sign
  std::swap(sign[i], sign[j]);

  // keep track of order of swaps
  swaps[i] = j;

  DataCollector::get()->countSwap();
}

void applySwaps(const int* swaps, int nrow, int ncol, double* R) {
  // apply the column swaps to block R
  for (int i = 0; i < ncol; ++i) {
    if (swaps[i] != i) {
      // swap col i and col swaps[i]
      callAndTime_dswap(nrow, &R[i], ncol, &R[swaps[i]], ncol);
    }
  }
}