#include "CholmodSolver.h"

void CholmodSolver::clear() {
  cholmod_free_factor(&this->L_, &this->c_);
  cholmod_free_sparse(&this->a_, &this->c_);
  cholmod_free_triplet(&this->T_, &this->c_);
  cholmod_free_dense(&this->x_, &this->c_);
  cholmod_free_dense(&this->b_, &this->c_);
  cholmod_finish(&c_);

  valid_ = false;
}

int CholmodSolver::factorAS(const HighsSparseMatrix& matrix,
                            const std::vector<double>& scaling) {
  std::cerr << "Cholmod does not factorize indefinite matrices\n";
  return kDecomposerStatusErrorFactorize;
}

int CholmodSolver::factorNE(const HighsSparseMatrix& highs_a,
                            const std::vector<double>& scaling) {
  // only execute factorization if it has not been done yet
  assert(!this->valid_);

  // compute normal equations matrix
  HighsSparseMatrix AThetaAT;
  int AAT_status = computeAThetaAT(highs_a, scaling, AThetaAT);
  if (AAT_status) return AAT_status;

  cholmod_start(&this->c_);

  // create the cholmod_triplet structure
  // cholmod_triplet* T = cholmod_allocate_triplet(AThetaAT.num_row_,
  // AThetaAT.num_col_, AThetaAT.value_.size(), 1, CHOLMOD_REAL,
  // &this->c);
  this->T_ = cholmod_allocate_triplet(AThetaAT.num_row_, AThetaAT.num_col_,
                                      AThetaAT.value_.size(), 1, CHOLMOD_REAL,
                                      &this->c_);
  int* Ti = static_cast<int*>(this->T_->i);
  int* Tj = static_cast<int*>(this->T_->j);
  double* Tx = static_cast<double*>(this->T_->x);
  // Extract upper triangular part of AThetaAT
  for (int col = 0; col < AThetaAT.num_col_; col++) {
    for (int idx = AThetaAT.start_[col]; idx < AThetaAT.start_[col + 1];
         idx++) {
      if (col >= AThetaAT.index_[idx]) {
        *Ti++ = AThetaAT.index_[idx];
        *Tj++ = col;
        *Tx++ = AThetaAT.value_[idx];
      }
    }
  }
  this->T_->nnz = Tx - static_cast<double*>(this->T_->x);

  //   cholmod_print_triplet(T, "T", &c);
  cholmod_check_triplet(this->T_, &this->c_);

  this->a_ = cholmod_triplet_to_sparse(this->T_, 1, &this->c_);

  //    cholmod_print_sparse(a,"A",&c);
  cholmod_check_sparse(this->a_, &this->c_);

  /*
  If you are going to factorize hundreds or more matrices with the same
  * nonzero pattern, you may wish to spend a great deal of time finding a
  * good permutation.  In this case, try setting Common->nmethods to 9.
  * The time spent in cholmod_analysis will be very high, but you need to
  * call it only once.
  */
  // c.nmethods = num_thods;
  // c.method[2].ordering = CHOLMOD_METIS;
  // c.postorder = 0;
  this->L_ = cholmod_analyze(this->a_, &this->c_);

  cholmod_factorize(this->a_, this->L_, &this->c_);

  size_t nnz = 0;
  if (this->L_->is_super) {  // Supernodal factor
    for (size_t s = 0; s < this->L_->nsuper; s++) {
      int pi = ((int*)(this->L_->pi))[s];
      int pj = ((int*)(this->L_->pi))[s + 1];
      for (int p = pi; p < pj; p++) {
        int i = ((int*)(this->L_->s))[p];
        nnz += (this->L_->n - i);
      }
    }
  } else {  // Simplicial factor
    for (size_t j = 0; j < this->L_->n; j++) {
      nnz += ((int*)(this->L_->nz))[j];
    }
  }

  this->valid_ = true;
  return kDecomposerStatusOk;
}

int CholmodSolver::solveNE(const std::vector<double>& rhs,
                           std::vector<double>& lhs) {
  // only execute the solve if factorization is valid
  assert(this->valid_);

  lhs = rhs;

  int system_size = lhs.size();

  // Create a cholmod_dense structure for rhs
  cholmod_dense* rhs_dense = cholmod_allocate_dense(system_size, 1, system_size,
                                                    CHOLMOD_REAL, &(this->c_));
  std::copy_n(lhs.data(), system_size, static_cast<double*>(rhs_dense->x));

  cholmod_dense* x = cholmod_solve(CHOLMOD_A, this->L_, rhs_dense, &this->c_);

  // Copy the solution back into rhs
  std::copy_n(static_cast<double*>(x->x), system_size, lhs.data());

  cholmod_free_dense(&rhs_dense, &this->c_);

  return kDecomposerStatusOk;
}

int CholmodSolver::solveAS(const std::vector<double>& rhs_x,
                           const std::vector<double>& rhs_y,
                           std::vector<double>& lhs_x,
                           std::vector<double>& lhs_y) {
  std::cerr << "Cholmod does not factorize indefinite matrices\n";
  return kDecomposerStatusErrorFactorize;
}