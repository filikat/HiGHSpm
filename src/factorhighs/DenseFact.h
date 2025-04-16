#ifndef FACTORHIGHS_DENSE_FACT_H
#define FACTORHIGHS_DENSE_FACT_H

/*
  Names:
  denseFact:
  - K : factorization kernel for diagonal blocks
  - F : blocked factorization in full format
  - FP: blocked factorization in format FP
  - FH : blocked factorization in "hybrid formats"

  Formats used:
  - F : Full format
  - P : lower Packed format
  - FP: lower Packed format with Full diagonal blocks
  - H : lower-blocked-Hybrid format
  - FH: lower-blocked-Hybrid format with Full diagonal blocks

  F, P do not use blocks. FP, H, FH use blocks.
  Blocks are always blocks of columns.
  F, P store by columns.
  FP stores by columns within the blocks. H, FH store by rows within the blocks.
  See report for more details.
*/

// dense factorization kernel
int denseFactK(char uplo, int n, double* A, int lda, int* pivot_sign,
               double thresh, double* regul, int* swaps, double* pivot_2x2,
               int sn, int bl, double max_in_R = -1);

// dense partial factorization, in full format
int denseFactF(int n, int k, int nb, double* A, int lda, double* B, int ldb,
               const int* pivot_sign, double thresh, double* regul, int sn);

// dense partial factorization, in packed format with full diagonal blocks
int denseFactFP(int n, int k, int nb, double* A, double* B,
                const int* pivot_sign, double thresh, double* regul, int sn);

// dense partial factorization, in "hybrid formats"
int denseFactFH(char format, int n, int k, int nb, double* A, double* B,
                const int* pivot_sign, double thresh, double* regul, int* swaps,
                double* pivot_2x2, int sn);

// function to convert A from lower packed, to lower-blocked-hybrid format
int denseFactFP2FH(double* A, int nrow, int ncol, int nb);

#endif