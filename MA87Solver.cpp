#include "MA87Solver.h"

void MA87Solver::Clear() {
  if (this->keep) {
    wrapper_ma87_finalise(&this->keep, &this->control);
  }
  valid = false;
}

int MA87Solver::FactorAS(const HighsSparseMatrix &matrix,
                         const std::vector<double> &theta) {

  std::cerr << "MA87 does not factorize indefinite matrices\n";
  return kDecomposerStatusErrorFactorize;
}

int MA87Solver::FactorNE(const HighsSparseMatrix &highs_a,
                         const std::vector<double> &theta) {

  // only execute factorization if it has not been done yet
  assert(!this->valid);

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
  this->order.resize(AThetaAT.num_row_);
  wrapper_mc68_default_control(&this->control_perm);
  int ord = 1;
  wrapper_mc68_order(ord, n, ptr.data(), row.data(), this->order.data(),
                     &this->control_perm, &this->info_perm);
  if (this->info_perm.flag < 0)
    return kDecomposerStatusErrorFactorize;

  // factorize with MA87
  wrapper_ma87_default_control(&this->control);
  wrapper_ma87_analyse(n, ptr.data(), row.data(), this->order.data(),
                       &this->keep, &this->control, &this->info);
  if (this->info.flag < 0)
    return kDecomposerStatusErrorFactorize;

  wrapper_ma87_factor(n, ptr.data(), row.data(), val.data(), this->order.data(),
                      &this->keep, &this->control, &this->info);
  if (this->info.flag < 0)
    return kDecomposerStatusErrorFactorize;

  this->valid = true;
  return kDecomposerStatusOk;
}

int MA87Solver::SolveNE(const HighsSparseMatrix &highs_a,
                        const std::vector<double> &theta,
                        const std::vector<double> &rhs,
                        std::vector<double> &lhs) {
  // only execute the solve if factorization is valid
  assert(this->valid);

  // initialize lhs with rhs
  lhs = rhs;
  int system_size = lhs.size();

  // solve with ma87
  wrapper_ma87_solve(0, 1, system_size, lhs.data(), this->order.data(),
                     &this->keep, &this->control, &this->info);
  if (this->info.flag < 0)
    return kDecomposerStatusErrorFactorize;

  return kDecomposerStatusOk;
}

int MA87Solver::SolveAS(const HighsSparseMatrix &highs_a,
                        const std::vector<double> &theta,
                        const std::vector<double> &rhs_x,
                        const std::vector<double> &rhs_y,
                        std::vector<double> &lhs_x,
                        std::vector<double> &lhs_y) {

  std::cerr << "MA87 does not factorize indefinite matrices\n";
  return kDecomposerStatusErrorFactorize;
}