#include "MA57Solver.h"

void MA57Solver::clear() {
  int info;
  wrapper_ma57_finalise(&this->factors_, &this->control_, &info);

  valid_ = false;
}

int MA57Solver::factorAS(const HighsSparseMatrix& matrix,
                         const std::vector<double>& theta) {
  std::cerr << "MA57 does not factorize indefinite matrices\n";
  return kDecomposerStatusErrorFactorize;
}

int MA57Solver::factorNE(const HighsSparseMatrix& highs_a,
                         const std::vector<double>& theta) {
  // only execute factorization if it has not been done yet
  assert(!this->valid_);

  // compute normal equations matrix
  HighsSparseMatrix AThetaAT;
  int AAT_status = computeAThetaAT(highs_a, theta, AThetaAT);
  if (AAT_status) return AAT_status;

  // initialize data to extract lower triangle
  int n = AThetaAT.num_col_;
  std::vector<int> ptr;

  // Extract lower triangular part of AAT
  for (int c = 0; c < n; c++) {
    ptr.push_back(val_.size());
    for (int idx = AThetaAT.start_[c]; idx < AThetaAT.start_[c + 1];
         idx++) {
      int row_idx = AThetaAT.index_[idx];
      if (row_idx >= c) {
        val_.push_back(AThetaAT.value_[idx]);
        row_.push_back(row_idx);
      }
    }
  }
  ptr.push_back(val_.size());

  int nz = ptr.back();

  // ma57 input is in triplet form
  col_.resize(ptr.back());
  for (int c = 0; c < n; ++c) {
    for (int el = ptr[c]; el < ptr[c + 1]; ++el) {
      col_[el] = c;
    }
  }

  wrapper_ma57_default_control(&this->control_);

  // switch off pivoting
  this->control_.pivoting = 3;
  this->control_.ordering = 2;

  wrapper_ma57_init_factors(&this->factors_);
  wrapper_ma57_analyse(n, nz, row_.data(), col_.data(), &this->factors_,
                       &this->control_, &this->ainfo_, this->order_.data());
  if (ainfo_.flag < 0) return kDecomposerStatusErrorFactorize;

  wrapper_ma57_factorize(n, nz, row_.data(), col_.data(), val_.data(),
                         &this->factors_, &this->control_, &this->finfo_);
  if (finfo_.flag < 0) return kDecomposerStatusErrorFactorize;

  this->valid_ = true;
  return kDecomposerStatusOk;
}

int MA57Solver::solveNE(const HighsSparseMatrix& highs_a,
                        const std::vector<double>& theta,
                        const std::vector<double>& rhs,
                        std::vector<double>& lhs) {
  // only execute the solve if factorization is valid
  assert(this->valid_);

  // initialize lhs with rhs
  lhs = rhs;
  int system_size = lhs.size();

  // solve with ma57
  wrapper_ma57_solve(system_size, row_.size(), row_.data(), col_.data(),
                     val_.data(), &this->factors_, lhs.data(), &this->control_,
                     &this->sinfo_);

  return kDecomposerStatusOk;
}

int MA57Solver::solveAS(const HighsSparseMatrix& highs_a,
                        const std::vector<double>& theta,
                        const std::vector<double>& rhs_x,
                        const std::vector<double>& rhs_y,
                        std::vector<double>& lhs_x,
                        std::vector<double>& lhs_y) {
  std::cerr << "MA57 does not factorize indefinite matrices\n";
  return kDecomposerStatusErrorFactorize;
}