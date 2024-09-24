#include "MA86Solver.h"

void MA86Solver::clear() {
  if (this->keep_) {
    wrapper_ma86_finalise(&this->keep_, &this->control_);
  }
  valid_ = false;
}

int MA86Solver::factorAS(const HighsSparseMatrix &matrix,
                         const std::vector<double> &theta) {

  // only execute factorization if it has not been done yet
  assert(!this->valid_);

  int n = matrix.num_col_;
  int m = matrix.num_row_;

  std::vector<int> ptr;
  std::vector<int> row;
  std::vector<double> val;

  // Extract lower triangular part of augmented system
  int row_index_offset = matrix.num_col_;
  for (int iCol = 0; iCol < matrix.num_col_; iCol++) {
    const double theta_i = !theta.empty() ? theta[iCol] : 1;

    ptr.push_back(val.size());
    val.push_back(-1 / theta_i);
    row.push_back(iCol);
    for (int iEl = matrix.start_[iCol]; iEl < matrix.start_[iCol + 1]; iEl++) {
      int iRow = matrix.index_[iEl];

      val.push_back(matrix.value_[iEl]);
      row.push_back(row_index_offset + iRow);
    }
  }

  const double diagonal = 0.0;
  for (int iRow = 0; iRow < matrix.num_row_; iRow++) {
    ptr.push_back(val.size());
    if (diagonal) {
      val.push_back(diagonal);
      row.push_back(row_index_offset + iRow);
    }
  }
  ptr.push_back(val.size());

  // ordering with MC68
  order_.resize(m + n);
  wrapper_mc68_default_control(&this->control_perm_);
  int ord = 1;
  wrapper_mc68_order(ord, n + m, ptr.data(), row.data(), this->order_.data(),
                     &this->control_perm_, &this->info_perm_);
  if (this->info_perm_.flag < 0)
    return kDecomposerStatusErrorFactorize;

  // factorization with MA86
  wrapper_ma86_default_control(&this->control_);
  wrapper_ma86_analyse(n + m, ptr.data(), row.data(), this->order_.data(),
                       &this->keep_, &this->control_, &this->info_);
  if (this->info_.flag < 0)
    return kDecomposerStatusErrorFactorize;

  wrapper_ma86_factor(n + m, ptr.data(), row.data(), val.data(),
                      this->order_.data(), &this->keep_, &this->control_,
                      &this->info_);
  if (this->info_.flag < 0)
    return kDecomposerStatusErrorFactorize;

  this->valid_ = true;
  return kDecomposerStatusOk;
}

int MA86Solver::factorNE(const HighsSparseMatrix &highs_a,
                         const std::vector<double> &theta) {

  // only execute factorization if it has not been done yet
  assert(!this->valid_);

  // compute normal equations matrix
  HighsSparseMatrix AThetaAT;
  int AAT_status = computeAThetaAT(highs_a, theta, AThetaAT);
  if (AAT_status)
    return AAT_status;

  // initialize data to extract lower triangle
  int n = AThetaAT.num_col_;
  const int array_base = 0;
  std::vector<int> ptr;
  std::vector<int> row;
  std::vector<double> val;

  // Extract lower triangular part of AAT
  for (int col = 0; col < AThetaAT.num_col_; col++) {
    ptr.push_back(val.size());
    for (int idx = AThetaAT.start_[col]; idx < AThetaAT.start_[col + 1];
         idx++) {
      int row_idx = AThetaAT.index_[idx];
      if (row_idx >= col) {
        val.push_back(AThetaAT.value_[idx]);
        row.push_back(row_idx);
      }
    }
  }
  ptr.push_back(val.size());

  // ordering with MC68
  this->order_.resize(AThetaAT.num_row_);
  wrapper_mc68_default_control(&this->control_perm_);
  int ord = 1;
  wrapper_mc68_order(ord, n, ptr.data(), row.data(), this->order_.data(),
                     &this->control_perm_, &this->info_perm_);
  if (this->info_perm_.flag < 0)
    return kDecomposerStatusErrorFactorize;

  // factorize with MA86
  wrapper_ma86_default_control(&this->control_);
  wrapper_ma86_analyse(n, ptr.data(), row.data(), this->order_.data(),
                       &this->keep_, &this->control_, &this->info_);
  if (this->info_.flag < 0)
    return kDecomposerStatusErrorFactorize;

  wrapper_ma86_factor(n, ptr.data(), row.data(), val.data(), this->order_.data(),
                      &this->keep_, &this->control_, &this->info_);
  if (this->info_.flag < 0)
    return kDecomposerStatusErrorFactorize;

  this->valid_ = true;
  return kDecomposerStatusOk;
}

int MA86Solver::solveNE(const HighsSparseMatrix &highs_a,
                        const std::vector<double> &theta,
                        const std::vector<double> &rhs,
                        std::vector<double> &lhs) {
  // only execute the solve if factorization is valid
  assert(this->valid_);

  // initialize lhs with rhs
  lhs = rhs;
  int system_size = lhs.size();

  // solve with ma86
  wrapper_ma86_solve(0, 1, system_size, lhs.data(), this->order_.data(),
                     &this->keep_, &this->control_, &this->info_);
  if (this->info_.flag < 0)
    return kDecomposerStatusErrorFactorize;

  return kDecomposerStatusOk;
}

int MA86Solver::solveAS(const HighsSparseMatrix &highs_a,
                        const std::vector<double> &theta,
                        const std::vector<double> &rhs_x,
                        const std::vector<double> &rhs_y,
                        std::vector<double> &lhs_x,
                        std::vector<double> &lhs_y) {

  // only execute the solve if factorization is valid
  assert(this->valid_);

  // create single rhs
  std::vector<double> rhs;
  rhs.insert(rhs.end(), rhs_x.begin(), rhs_x.end());
  rhs.insert(rhs.end(), rhs_y.begin(), rhs_y.end());

  int system_size = highs_a.num_col_ + highs_a.num_row_;

  // solve using ma86
  wrapper_ma86_solve(0, 1, system_size, rhs.data(), this->order_.data(),
                     &this->keep_, &this->control_, &this->info_);
  if (this->info_.flag < 0)
    return kDecomposerStatusErrorFactorize;

  // split lhs
  lhs_x = std::vector<double>(rhs.begin(), rhs.begin() + highs_a.num_col_);
  lhs_y = std::vector<double>(rhs.begin() + highs_a.num_col_, rhs.end());

  return kDecomposerStatusOk;
}