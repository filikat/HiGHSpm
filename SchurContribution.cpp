#include "Metis_caller.h"

// Functions to compute contributions to the schur complement

// Solve rhs one at a time.
// I keep this only for comparison, because it is slow
void Metis_caller::SchurContributionSingle(int i) {
  // to build Schur complement:
  // for each linking block
  // - access rows of linking block B^T (stored row-wise)
  // - create dense vector with row b
  // - do a forward solve and store as sparse column
  //    L^-1 * b
  // - do a diagonal solve and store as sparse column
  //    D^-1 * L^-1 * b
  // - leave L^-1 * B^T col-wise, but get D^-1 * L^-1 * B^T row-wise
  // - use same trick as computeAThetaAT to find contribution to Schur
  //    complement
  // - Store everything in a dense matrix

  int linkSize = blockSize.back();

  // current linking block
  HighsSparseMatrix& B = Blocks[2 * i + 1];

  // row-wise linking block
  HighsSparseMatrix Bt = B;
  Bt.ensureRowwise();

  // space for L^-1 B^T and D^-1 L^-1 B^T
  HighsSparseMatrix forwardSolvedBlock;
  HighsSparseMatrix diagForwardSolvedBlock;
  forwardSolvedBlock.num_row_ = B.num_col_;
  forwardSolvedBlock.num_col_ = B.num_row_;
  diagForwardSolvedBlock.num_row_ = B.num_col_;
  diagForwardSolvedBlock.num_col_ = B.num_row_;

  // go through the rows of B
  for (int row = 0; row < B.num_row_; ++row) {
    // space for dense vector
    std::vector<double> denseRow(B.num_col_, 0.0);

    // avoid computations if row is empty
    if (Bt.start_[row] == Bt.start_[row + 1]) {
      forwardSolvedBlock.start_.push_back(forwardSolvedBlock.start_.back());
      diagForwardSolvedBlock.start_.push_back(
          diagForwardSolvedBlock.start_.back());
      continue;
    }

    // go through the entries of the row and fill in dense_row
    for (int elem = Bt.start_[row]; elem < Bt.start_[row + 1]; ++elem) {
      denseRow[Bt.index_[elem]] = Bt.value_[elem];
    }

    double t1 = getWallTime();
    // on return from diagonalForwardSolve, denseRow contains L^-1 * denseRow
    // and diagForwardSolvedRow contains D^-1 * L^-1 * denseRow
    std::vector<double> diagForwardSolvedRow(B.num_col_, 0.0);
    diagonalForwardSolve(denseRow.data(), 1, invertData[i], expData[i],
                         diagForwardSolvedRow.data());
    schur_dfsolve_time += getWallTime() - t1;

    t1 = getWallTime();
    // transform denseRow into sparse column of forwardSolvedBlock
    int countNnzRow{};
    for (int i = 0; i < denseRow.size(); ++i) {
      if (std::fabs(denseRow[i]) > 1e-20) {
        forwardSolvedBlock.index_.push_back(i);
        forwardSolvedBlock.value_.push_back(denseRow[i]);
        ++countNnzRow;
      }
    }
    forwardSolvedBlock.start_.push_back(forwardSolvedBlock.start_.back() +
                                        countNnzRow);

    // transform diagForwardSolvedRow into sparse column of
    // diagForwardSolvedBlock
    countNnzRow = 0;
    for (int i = 0; i < diagForwardSolvedRow.size(); ++i) {
      if (std::fabs(diagForwardSolvedRow[i]) > 1e-20) {
        diagForwardSolvedBlock.index_.push_back(i);
        diagForwardSolvedBlock.value_.push_back(diagForwardSolvedRow[i]);
        ++countNnzRow;
      }
    }
    diagForwardSolvedBlock.start_.push_back(
        diagForwardSolvedBlock.start_.back() + countNnzRow);

    schur_transform_time += getWallTime() - t1;
  }

  // diagForwardSolvedBlock needs row-wise access
  diagForwardSolvedBlock.ensureRowwise();

  // similar trick as computeAThetaAT, with simplifications, to compute
  //      Q^T        *          R
  // (L^-1 * B^T)^T  *  (D^-1 * L^-1 * B^T)
  // ^^^^^^^^^^^^       ^^^^^^^^^^^^^^^^^^^
  //  col-wise               row-wise
  // access to Q            access to R
  //
  // Q and R may have different sparsity pattern, because D may have 2x2
  // pivots. This may cause inefficiencies.

  HighsSparseMatrix& Q = forwardSolvedBlock;
  HighsSparseMatrix& R = diagForwardSolvedBlock;

  double t1 = getWallTime();

  // go through columns of Q
  for (int col = 0; col < Q.num_col_; ++col) {
    // go through entries of the column
    for (int colEl = Q.start_[col]; colEl < Q.start_[col + 1]; ++colEl) {
      // row index of current entry
      int row = Q.index_[colEl];
      // go through entries of the row of R
      for (int rowEl = R.start_[row]; rowEl < R.start_[row + 1]; ++rowEl) {
        // this entry contributes to Schur(col,col1) and Schur (col1,col)
        int col1 = R.index_[rowEl];

        // avoid counting entries twice
        if (col1 < col) continue;

        double value = Q.value_[colEl] * R.value_[rowEl];
        lapack_a[col * linkSize + col1] -= value;
        if (col != col1) {
          lapack_a[col1 * linkSize + col] -= value;
        }
      }
    }
  }
  schur_multiply_time += getWallTime() - t1;
}

// Solve rhs all at once.
void Metis_caller::SchurContributionMultiple(int i) {
  // to build Schur complement:
  // for each linking block
  // - access rows of linking block B^T (stored row-wise)
  // - create dense vector with row b
  // - do a forward solve and store as sparse column
  //    L^-1 * b
  // - do a diagonal solve and store as sparse column
  //    D^-1 * L^-1 * b
  // - do the previous solve operations all at once
  // - leave L^-1 * B^T col-wise, but get D^-1 * L^-1 * B^T row-wise
  // - use same trick as computeAThetaAT to find contribution to Schur
  //    complement
  // - Store everything in a dense matrix
  int linkSize = blockSize.back();

  // current linking block
  HighsSparseMatrix& B = Blocks[2 * i + 1];

  // row-wise linking block
  HighsSparseMatrix Bt = B;
  Bt.ensureRowwise();

  // space for L^-1 B^T and D^-1 L^-1 B^T
  HighsSparseMatrix forwardSolvedBlock;
  HighsSparseMatrix diagForwardSolvedBlock;
  forwardSolvedBlock.num_row_ = B.num_col_;
  forwardSolvedBlock.num_col_ = B.num_row_;
  diagForwardSolvedBlock.num_row_ = B.num_col_;
  diagForwardSolvedBlock.num_col_ = B.num_row_;

  // space for all dense rows together
  std::vector<double> denseRows(B.num_col_ * B.num_row_, 0.0);

  // go through the rows of B and fill denseRows
  for (int row = 0; row < B.num_row_; ++row) {
    // offset for writing to denseRows
    int offset = B.num_col_ * row;

    // go through the entries of the row and fill in dense_rows
    for (int elem = Bt.start_[row]; elem < Bt.start_[row + 1]; ++elem) {
      denseRows[offset + Bt.index_[elem]] = Bt.value_[elem];
    }
  }

  double t1 = getWallTime();
  // diagonalForwardSolve of all rows at the same time.
  // On return from diagonalForwardSolve, denseRows contains L^-1 * denseRows
  // and diagForwardSolvedRows contains D^-1 * L^-1 * denseRows
  std::vector<double> diagForwardSolvedRows(B.num_col_ * B.num_row_, 0.0);
  if (linkSize > 0) {
    diagonalForwardSolve(denseRows.data(), B.num_row_, invertData[i],
                         expData[i], diagForwardSolvedRows.data());
  }
  schur_dfsolve_time += getWallTime() - t1;

  // go through the rows of B
  for (int row = 0; row < B.num_row_; ++row) {
    // offset for reading from denseRows
    int offset = B.num_col_ * row;

    t1 = getWallTime();
    // transform denseRows into sparse columns of forwardSolvedBlock
    int countNnzRow{};
    for (int i = 0; i < B.num_col_; ++i) {
      if (std::fabs(denseRows[offset + i]) > 1e-20) {
        forwardSolvedBlock.index_.push_back(i);
        forwardSolvedBlock.value_.push_back(denseRows[offset + i]);
        ++countNnzRow;
      }
    }
    forwardSolvedBlock.start_.push_back(forwardSolvedBlock.start_.back() +
                                        countNnzRow);

    // transform diagForwardSolvedRows into sparse columns of
    // diagForwardSolvedBlock
    countNnzRow = 0;
    for (int i = 0; i < B.num_col_; ++i) {
      if (std::fabs(diagForwardSolvedRows[offset + i]) > 1e-20) {
        diagForwardSolvedBlock.index_.push_back(i);
        diagForwardSolvedBlock.value_.push_back(
            diagForwardSolvedRows[offset + i]);
        ++countNnzRow;
      }
    }
    diagForwardSolvedBlock.start_.push_back(
        diagForwardSolvedBlock.start_.back() + countNnzRow);

    schur_transform_time += getWallTime() - t1;
  }

  // diagForwardSolvedBlock needs row-wise access
  diagForwardSolvedBlock.ensureRowwise();

  // similar trick as computeAThetaAT, with simplifications, to compute
  //      Q^T        *          R
  // (L^-1 * B^T)^T  *  (D^-1 * L^-1 * B^T)
  // ^^^^^^^^^^^^       ^^^^^^^^^^^^^^^^^^^
  //  col-wise               row-wise
  // access to Q            access to R
  //
  // Q and R may have different sparsity pattern, because D may have 2x2
  // pivots. This may cause inefficiencies.

  HighsSparseMatrix& Q = forwardSolvedBlock;
  HighsSparseMatrix& R = diagForwardSolvedBlock;

  t1 = getWallTime();

  // go through columns of Q
  for (int col = 0; col < Q.num_col_; ++col) {
    // go through entries of the column
    for (int colEl = Q.start_[col]; colEl < Q.start_[col + 1]; ++colEl) {
      // row index of current entry
      int row = Q.index_[colEl];
      // go through entries of the row of R
      for (int rowEl = R.start_[row]; rowEl < R.start_[row + 1]; ++rowEl) {
        // this entry contributes to Schur(col,col1) and Schur (col1,col)
        int col1 = R.index_[rowEl];

        // avoid counting entries twice
        if (col1 < col) continue;

        double value = Q.value_[colEl] * R.value_[rowEl];

        // insert element in position (col, col1)
        lapack_a[col * linkSize + col1] -= value;
        if (col != col1) {
          // insert element in position (col1, col)
          lapack_a[col1 * linkSize + col] -= value;
        }
      }
    }
  }
  schur_multiply_time += getWallTime() - t1;
}

// Use Highs LU factorization
void Metis_caller::SchurContributionHFactor(int part) {}
