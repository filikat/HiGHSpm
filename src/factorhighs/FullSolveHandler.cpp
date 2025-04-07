#include "FullSolveHandler.h"
#include "Auxiliary.h"
#include "CallAndTimeBlas.h"
#include "DataCollector.h"

FullSolveHandler::FullSolveHandler(
    const Symbolic& S, const std::vector<std::vector<double>>& sn_columns)
    : SolveHandler(S, sn_columns) {}

void FullSolveHandler::forwardSolve(std::vector<double>& x) const {
  // Forward solve.
  // Blas calls: dtrsv, dgemv

  // supernode columns in format F

  for (int sn = 0; sn < S_.sn(); ++sn) {
    // leading size of supernode
    const int ldSn = S_.ptr(sn + 1) - S_.ptr(sn);

    // number of columns in the supernode
    const int sn_size = S_.snStart(sn + 1) - S_.snStart(sn);

    // first colums of the supernode
    const int sn_start = S_.snStart(sn);

    // size of clique of supernode
    const int clique_size = ldSn - sn_size;

    // index to access S->rows for this supernode
    const int start_row = S_.ptr(sn);

    callAndTime_dtrsv('L', 'N', 'U', sn_size, sn_columns_[sn].data(), ldSn,
                      &x[sn_start], 1);

    // temporary space for gemv
    std::vector<double> y(clique_size);

    callAndTime_dgemv('N', clique_size, sn_size, 1.0, &sn_columns_[sn][sn_size],
                      ldSn, &x[sn_start], 1, 0.0, y.data(), 1);

    // scatter solution of gemv
    for (int i = 0; i < clique_size; ++i) {
      const int row = S_.rows(start_row + sn_size + i);
      x[row] -= y[i];
    }
  }
}

void FullSolveHandler::backwardSolve(std::vector<double>& x) const {
  // Backward solve.
  // Blas calls: dtrsv, dgemv

  // supernode columns in format F

  // go through the sn in reverse order
  for (int sn = S_.sn() - 1; sn >= 0; --sn) {
    // leading size of supernode
    const int ldSn = S_.ptr(sn + 1) - S_.ptr(sn);

    // number of columns in the supernode
    const int sn_size = S_.snStart(sn + 1) - S_.snStart(sn);

    // first colums of the supernode
    const int sn_start = S_.snStart(sn);

    // size of clique of supernode
    const int clique_size = ldSn - sn_size;

    // index to access S->rows for this supernode
    const int start_row = S_.ptr(sn);

    // temporary space for gemv
    std::vector<double> y(clique_size);

    // scatter entries into y
    for (int i = 0; i < clique_size; ++i) {
      const int row = S_.rows(start_row + sn_size + i);
      y[i] = x[row];
    }

    callAndTime_dgemv('T', clique_size, sn_size, -1.0,
                      &sn_columns_[sn][sn_size], ldSn, y.data(), 1, 1.0,
                      &x[sn_start], 1);

    callAndTime_dtrsv('L', 'T', 'U', sn_size, sn_columns_[sn].data(), ldSn,
                      &x[sn_start], 1);
  }
}

void FullSolveHandler::diagSolve(std::vector<double>& x) const {
  // Diagonal solve

  // supernode columns in format F

  for (int sn = 0; sn < S_.sn(); ++sn) {
    // leading size of supernode
    const int ldSn = S_.ptr(sn + 1) - S_.ptr(sn);

    for (int col = S_.snStart(sn); col < S_.snStart(sn + 1); ++col) {
      // relative index of column within supernode
      const int j = col - S_.snStart(sn);

      // diagonal entry of column j
      const double d = sn_columns_[sn][j + j * ldSn];

      x[col] /= d;
    }
  }
}