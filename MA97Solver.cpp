#include "MA97Solver.h"

void MA97Solver::clear() {
  if (this->akeep_ && this->fkeep_) {
    wrapper_ma97_finalise(&this->akeep_, &this->fkeep_);
  }
  valid_ = false;
}

int MA97Solver::factorAS(const HighsSparseMatrix& matrix,
                         const std::vector<double>& scaling) {
  std::cerr << "MA97 does not factorize indefinite matrices\n";
  return kDecomposerStatusErrorFactorize;
}

int MA97Solver::factorNE(const HighsSparseMatrix& highs_a,
                         const std::vector<double>& scaling) {
  // only execute factorization if it has not been done yet
  assert(!this->valid_);

  // compute normal equations matrix
  HighsSparseMatrix AThetaAT;
  int AAT_status = computeAThetaAT(highs_a, scaling, AThetaAT);
  if (AAT_status) return AAT_status;

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
  if (this->info_perm_.flag < 0) return kDecomposerStatusErrorFactorize;

  // factorize with MA97
  wrapper_ma97_default_control(&this->control_);
  wrapper_ma97_analyse(n, ptr.data(), row.data(), this->order_.data(),
                       &this->akeep_, &this->control_, &this->info_);
  if (this->info_.flag < 0) return kDecomposerStatusErrorFactorize;

  wrapper_ma97_factor(3, n, ptr.data(), row.data(), val.data(), &this->akeep_,
                      &this->fkeep_, &this->control_, &this->info_);
  if (this->info_.flag < 0) return kDecomposerStatusErrorFactorize;

  this->valid_ = true;
  return kDecomposerStatusOk;
}

int MA97Solver::solveNE(const std::vector<double>& rhs,
                        std::vector<double>& lhs) {
  // only execute the solve if factorization is valid
  assert(this->valid_);

  // initialize lhs with rhs
  lhs = rhs;
  int system_size = lhs.size();

  // solve with ma97
  wrapper_ma97_solve(0, 1, system_size, lhs.data(), &this->akeep_,
                     &this->fkeep_, &this->control_, &this->info_);
  if (this->info_.flag < 0) return kDecomposerStatusErrorFactorize;

  return kDecomposerStatusOk;
}

int MA97Solver::solveAS(const std::vector<double>& rhs_x,
                        const std::vector<double>& rhs_y,
                        std::vector<double>& lhs_x,
                        std::vector<double>& lhs_y) {
  std::cerr << "MA97 does not factorize indefinite matrices\n";
  return kDecomposerStatusErrorFactorize;
}