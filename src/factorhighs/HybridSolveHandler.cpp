#include "HybridSolveHandler.h"

#include "FactorHiGHSSettings.h"
#include "FormatHandler.h"
#include "Auxiliary.h"
#include "CallAndTimeBlas.h"
#include "DataCollector.h"

HybridSolveHandler::HybridSolveHandler(
    const Symbolic& S, const std::vector<std::vector<double>>& sn_columns,
    const std::vector<std::vector<int>>& swaps,
    const std::vector<std::vector<double>>& pivot_2x2)
    : SolveHandler(S, sn_columns), swaps_{swaps}, pivot_2x2_{pivot_2x2} {}

void HybridSolveHandler::forwardSolve(std::vector<double>& x) const {
  // Forward solve.
  // Blas calls: dtrsv, dgemv

  // supernode columns in format FH

  const int nb = S_.blockSize();

  for (int sn = 0; sn < S_.sn(); ++sn) {
    // leading size of supernode
    const int ldSn = S_.ptr(sn + 1) - S_.ptr(sn);

    // number of columns in the supernode
    const int sn_size = S_.snStart(sn + 1) - S_.snStart(sn);

    // first colums of the supernode
    const int sn_start = S_.snStart(sn);

    // index to access S->rows for this supernode
    const int start_row = S_.ptr(sn);

    // number of blocks of columns
    const int n_blocks = (sn_size - 1) / nb + 1;

    // index to access snColumns[sn]
    int SnCol_ind{};

    // go through blocks of columns for this supernode
    for (int j = 0; j < n_blocks; ++j) {
      // number of columns in the block
      const int jb = std::min(nb, sn_size - nb * j);

      // number of entries in diagonal part
      const int diag_entries = jb * jb;

      // index to access vector x
      const int x_start = sn_start + nb * j;

#ifdef PIVOTING
      // apply swaps to portion of rhs that is affected
      const int* current_swaps = &swaps_[sn][nb * j];
      permuteWithSwaps(&x[x_start], current_swaps, jb);
#endif

      callAndTime_dtrsv('U', 'T', 'U', jb, &sn_columns_[sn][SnCol_ind], jb,
                        &x[x_start], 1);

      SnCol_ind += diag_entries;

      // temporary space for gemv
      const int gemv_space = ldSn - nb * j - jb;
      std::vector<double> y(gemv_space);

      callAndTime_dgemv('T', jb, gemv_space, 1.0, &sn_columns_[sn][SnCol_ind],
                        jb, &x[x_start], 1, 0.0, y.data(), 1);
      SnCol_ind += jb * gemv_space;

      // scatter solution of gemv
      for (int i = 0; i < gemv_space; ++i) {
        const int row = S_.rows(start_row + nb * j + jb + i);
        x[row] -= y[i];
      }

#ifdef PIVOTING
      // apply inverse swaps
      permuteWithSwaps(&x[x_start], current_swaps, jb, true);
#endif
    }
  }
}

void HybridSolveHandler::backwardSolve(std::vector<double>& x) const {
  // Backward solve.
  // Blas calls: dtrsv, dgemv

  // supernode columns in format FH

  const int nb = S_.blockSize();

  // go through the sn in reverse order
  for (int sn = S_.sn() - 1; sn >= 0; --sn) {
    // leading size of supernode
    const int ldSn = S_.ptr(sn + 1) - S_.ptr(sn);

    // number of columns in the supernode
    const int sn_size = S_.snStart(sn + 1) - S_.snStart(sn);

    // first colums of the supernode
    const int sn_start = S_.snStart(sn);

    // index to access S->rows for this supernode
    const int start_row = S_.ptr(sn);

    // number of blocks of columns
    const int n_blocks = (sn_size - 1) / nb + 1;

    // index to access snColumns[sn]
    // initialized with the total number of entries of snColumns[sn]
    int SnCol_ind = sn_columns_[sn].size() - extra_space;

    // go through blocks of columns for this supernode in reverse order
    for (int j = n_blocks - 1; j >= 0; --j) {
      // number of columns in the block
      const int jb = std::min(nb, sn_size - nb * j);

      // number of entries in diagonal part
      const int diag_entries = jb * jb;

      // index to access vector x
      const int x_start = sn_start + nb * j;

#ifdef PIVOTING
      // apply swaps to portion of rhs that is affected
      const int* current_swaps = &swaps_[sn][nb * j];
      permuteWithSwaps(&x[x_start], current_swaps, jb);
#endif

      // temporary space for gemv
      const int gemv_space = ldSn - nb * j - jb;
      std::vector<double> y(gemv_space);

      // scatter entries into y
      for (int i = 0; i < gemv_space; ++i) {
        const int row = S_.rows(start_row + nb * j + jb + i);
        y[i] = x[row];
      }

      SnCol_ind -= jb * gemv_space;
      callAndTime_dgemv('N', jb, gemv_space, -1.0, &sn_columns_[sn][SnCol_ind],
                        jb, y.data(), 1, 1.0, &x[x_start], 1);

      SnCol_ind -= diag_entries;
      callAndTime_dtrsv('U', 'N', 'U', jb, &sn_columns_[sn][SnCol_ind], jb,
                        &x[x_start], 1);

#ifdef PIVOTING
      // apply inverse swaps
      permuteWithSwaps(&x[x_start], current_swaps, jb, true);
#endif
    }
  }
}

void HybridSolveHandler::diagSolve(std::vector<double>& x) const {
  // Diagonal solve

  // supernode columns in format FH

  const int nb = S_.blockSize();

  for (int sn = 0; sn < S_.sn(); ++sn) {
    // leading size of supernode
    const int ldSn = S_.ptr(sn + 1) - S_.ptr(sn);

    // number of columns in the supernode
    const int sn_size = S_.snStart(sn + 1) - S_.snStart(sn);

    // first colums of the supernode
    const int sn_start = S_.snStart(sn);

    // number of blocks of columns
    const int n_blocks = (sn_size - 1) / nb + 1;

    // index to access diagonal part of block
    int diag_start{};

    // go through blocks of columns for this supernode
    for (int j = 0; j < n_blocks; ++j) {
      // number of columns in the block
      const int jb = std::min(nb, sn_size - nb * j);

#ifdef PIVOTING
      // apply swaps to portion of rhs that is affected
      const int* current_swaps = &swaps_[sn][nb * j];
      permuteWithSwaps(&x[sn_start + nb * j], current_swaps, jb);
#endif

      const double* current_2x2 = &pivot_2x2_[sn][nb * j];
      int step = 1;

      // go through columns of block
      for (int col = 0; col < jb; col += step) {
        if (current_2x2[col] == 0.0) {
          // 1x1 pivots
          step = 1;
          const double inv_d = sn_columns_[sn][diag_start + col + jb * col];
          x[sn_start + nb * j + col] *= inv_d;
        } else {
          // 2x2 pivots
          step = 2;

          // inverse of 2x2 pivot
          const double i_d1 = sn_columns_[sn][diag_start + col + jb * col];
          const double i_d2 =
              sn_columns_[sn][diag_start + col + 1 + jb * (col + 1)];
          const double i_off = current_2x2[col];

          double x1 = x[sn_start + nb * j + col];
          double x2 = x[sn_start + nb * j + col + 1];

          x[sn_start + nb * j + col] = i_d1 * x1 + i_off * x2;
          x[sn_start + nb * j + col + 1] = i_d2 * x2 + i_off * x1;
        }
      }

#ifdef PIVOTING
      // apply inverse swaps
      permuteWithSwaps(&x[sn_start + nb * j], current_swaps, jb, true);
#endif

      // move diag_start forward by number of diagonal entries in block
      diag_start += jb * jb;

      // move diag_start forward by number of sub-diagonal entries in block
      diag_start += (ldSn - nb * j - jb) * jb;
    }
  }
}