#include "HFactorSolver.h"

void HFactorSolver::clear() {
  this->basic_index_.clear();
  valid_ = false;
}

int HFactorSolver::factorAS(const HighsSparseMatrix& matrix,
                            const std::vector<double>& scaling) {
  // only execute factorization if it has not been done yet
  assert(!this->valid_);

  // Prepare data structures for HiGHS
  HighsSparseMatrix augmented;
  // Extract lower triangular part of augmented system
  //
  // First the columns for [-\Theta^{-1}]
  //                       [      A     ]
  int row_index_offset = matrix.num_col_;
  std::vector<HighsInt>& start = augmented.start_;
  std::vector<HighsInt>& index = augmented.index_;
  std::vector<double>& value = augmented.value_;
  std::vector<HighsInt>& basic_index = this->basic_index_;
  std::vector<double> max_scaled_matrix_row_value(matrix.num_row_, 0);
  HFactor& factor = this->factor_;
  std::vector<HighsInt> AT_start(matrix.num_row_, 0);
  assert(basic_index.size() == 0);
  for (int iCol = 0; iCol < matrix.num_col_; iCol++) {
    basic_index.push_back(iCol);
    const double scaling_i = !scaling.empty() ? scaling[iCol] : 1;
    double diag = -scaling_i;
    assert(diag < 0);
    value.push_back(diag);
    index.push_back(iCol);
    for (int iEl = matrix.start_[iCol]; iEl < matrix.start_[iCol + 1]; iEl++) {
      const int iRow = matrix.index_[iEl];
      AT_start[iRow]++;
      value.push_back(matrix.value_[iEl]);
      index.push_back(row_index_offset + iRow);
    }
    start.push_back(index.size());
  }
  // Now the columns for [A^T]
  //                     [ 0 ]
  //
  // Set the starts for the remaining columns, and initialise AT_start
  // for inserting indices and values
  for (int iRow = 0; iRow < matrix.num_row_; iRow++) {
    basic_index.push_back(row_index_offset + iRow);
    start.push_back(start[row_index_offset + iRow] + AT_start[iRow]);
    AT_start[iRow] = start[row_index_offset + iRow];
  }
  const int system_size = matrix.num_col_ + matrix.num_row_;
  const int augmented_nnz = start[system_size];
  index.resize(augmented_nnz);
  value.resize(augmented_nnz);
  // Work through the columns of A as rows of A^T
  for (int iCol = 0; iCol < matrix.num_col_; iCol++) {
    for (int iEl = matrix.start_[iCol]; iEl < matrix.start_[iCol + 1]; iEl++) {
      const int iRow = matrix.index_[iEl];
      index[AT_start[iRow]] = iCol;
      value[AT_start[iRow]] = matrix.value_[iEl];
      AT_start[iRow]++;
    }
  }
  augmented.num_col_ = system_size;
  augmented.num_row_ = system_size;
  assert(int(basic_index.size()) == system_size);

  factor.setup(augmented, basic_index);

  const HighsInt rank_deficiency = factor.build();

  int invert_status;
  if (rank_deficiency) {
    invert_status = kDecomposerStatusErrorFactorize;
  } else {
    this->valid_ = true;
    invert_status = kDecomposerStatusOk;
  }
  return invert_status;
}

int HFactorSolver::factorNE(const HighsSparseMatrix& highs_a,
                            const std::vector<double>& scaling) {
  // only execute factorization if it has not been done yet
  assert(!this->valid_);

  // compute normal equations matrix
  HighsSparseMatrix AThetaAT;
  int AAT_status = computeAThetaAT(highs_a, scaling, AThetaAT);
  if (AAT_status) return AAT_status;

  HFactor& factor = this->factor_;
  std::vector<HighsInt>& basic_index = this->basic_index_;
  for (int iCol = 0; iCol < AThetaAT.num_col_; iCol++)
    basic_index.push_back(iCol);

  factor.setup(AThetaAT, basic_index);

  const HighsInt rank_deficiency = factor.build();

  int invert_status;
  if (rank_deficiency) {
    invert_status = kDecomposerStatusErrorFactorize;
  } else {
    this->valid_ = true;
    invert_status = kDecomposerStatusOk;
  }
  return invert_status;
}

int HFactorSolver::solveNE(const std::vector<double>& rhs,
                           std::vector<double>& lhs) {
  // only execute the solve if factorization is valid
  assert(this->valid_);

  std::vector<double> sol = rhs;
  this->factor_.ftranCall(sol);
  for (int iCol = 0; iCol < int(rhs.size()); iCol++)
    lhs[this->basic_index_[iCol]] = sol[iCol];

  return kDecomposerStatusOk;
}

int HFactorSolver::solveAS(const std::vector<double>& rhs_x,
                           const std::vector<double>& rhs_y,
                           std::vector<double>& lhs_x,
                           std::vector<double>& lhs_y) {
  // only execute the solve if factorization is valid
  assert(this->valid_);

  // create single rhs
  std::vector<double> rhs;
  rhs.insert(rhs.end(), rhs_x.begin(), rhs_x.end());
  rhs.insert(rhs.end(), rhs_y.begin(), rhs_y.end());

  int m = rhs_y.size();
  int n = rhs_x.size();
  int system_size = m + n;

  // solve using HFactor
  std::vector<double> sol = rhs;
  this->factor_.ftranCall(sol);
  for (int iCol = 0; iCol < int(rhs.size()); iCol++)
    rhs[this->basic_index_[iCol]] = sol[iCol];

  // split lhs
  lhs_x = std::vector<double>(rhs.begin(), rhs.begin() + n);
  lhs_y = std::vector<double>(rhs.begin() + n, rhs.end());

  return kDecomposerStatusOk;
}