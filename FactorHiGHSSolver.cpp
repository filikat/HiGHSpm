#include "FactorHiGHSSolver.h"

void FactorHiGHSSolver::clear() {
  valid_ = false;
  worst_res_ = 0.0;
}

int FactorHiGHSSolver::setup(const HighsSparseMatrix& A,
                             const std::vector<int>& parameters) {
  std::vector<int> ptrLower;
  std::vector<int> rowsLower;

  int nA = A.num_col_;
  int mA = A.num_row_;
  int nzA = A.numNz();

  int negative_pivots{};

  S_.setFact((FactType)parameters[kParamFact]);
  S_.setFormat((FormatType)parameters[kParamFormat]);

  int nla_type = parameters[kParamNla];

  // Build the matrix
  if (nla_type == kOptionNlaAugmented) {
    // Augmented system, lower triangular

    ptrLower.resize(nA + mA + 1);
    rowsLower.resize(nA + nzA + mA);

    int next = 0;

    for (int i = 0; i < nA; ++i) {
      // diagonal element
      rowsLower[next] = i;
      ++next;

      // column of A
      for (int el = A.start_[i]; el < A.start_[i + 1]; ++el) {
        rowsLower[next] = A.index_[el] + nA;
        ++next;
      }

      ptrLower[i + 1] = next;
    }

    // 2,2 block
    for (int i = 0; i < mA; ++i) {
      rowsLower[next] = nA + i;
      ++next;
      ptrLower[nA + i + 1] = ptrLower[nA + i] + 1;
    }

    negative_pivots = nA;

  } else {
    // Normal equations, full matrix
    std::vector<double> theta;
    HighsSparseMatrix AAt;
    int status = computeLowerAThetaAT(A, theta, AAt);
    if (status) {
      printf("Failure: AAt is too large\n");
      return kRetOutOfMemory;
    }

    rowsLower = std::move(AAt.index_);
    ptrLower = std::move(AAt.start_);
  }

  // Perform analyse phase
  Analyse analyse(rowsLower, ptrLower, {}, negative_pivots);
  analyse.run(S_);
  // S_.print();

  return kRetOk;
}

int FactorHiGHSSolver::factorAS(const HighsSparseMatrix& A,
                                const std::vector<double>& theta) {
  // only execute factorization if it has not been done yet
  assert(!this->valid_);

  // initialize
  std::vector<int> ptrLower;
  std::vector<int> rowsLower;
  std::vector<double> valLower;

  int nA = A.num_col_;
  int mA = A.num_row_;
  int nzA = A.numNz();

  ptrLower.resize(nA + mA + 1);
  rowsLower.resize(nA + nzA + mA);
  valLower.resize(nA + nzA + mA);

  // build lower triangle
  int next = 0;

  for (int i = 0; i < nA; ++i) {
    // diagonal element
    rowsLower[next] = i;
    valLower[next++] = -1.0 / theta[i];

    // column of A
    for (int el = A.start_[i]; el < A.start_[i + 1]; ++el) {
      rowsLower[next] = A.index_[el] + nA;
      valLower[next++] = A.value_[el];
    }

    ptrLower[i + 1] = next;
  }

  // 2,2 block
  for (int i = 0; i < mA; ++i) {
    rowsLower[next] = nA + i;
    valLower[next++] = kDualStaticRegularization;
    ptrLower[nA + i + 1] = ptrLower[nA + i] + 1;
  }

  // factorise matrix
  Factorise factorise(S_, rowsLower, ptrLower, valLower);
  int status = factorise.run(N_);
  if (status) {
    return kDecomposerStatusErrorFactorize;
  }

  // save extreme values of the factorisation
  minD_ = factorise.minD_;
  maxD_ = factorise.maxD_;
  minoffD_ = factorise.minoffD_;
  maxoffD_ = factorise.maxoffD_;

  this->valid_ = true;
  use_as_ = true;
  return kDecomposerStatusOk;
}

int FactorHiGHSSolver::factorNE(const HighsSparseMatrix& A,
                                const std::vector<double>& theta) {
  // only execute factorization if it has not been done yet
  assert(!this->valid_);

  // initialize
  std::vector<int> ptrLower;
  std::vector<int> rowsLower;
  std::vector<double> valLower;

  int nA = A.num_col_;
  int mA = A.num_row_;
  int nzA = A.numNz();

  // build full matrix
  HighsSparseMatrix AAt;
  int status = computeLowerAThetaAT(A, theta, AAt);

  // factorise
  Factorise factorise(S_, AAt.index_, AAt.start_, AAt.value_);
  int factorise_status = factorise.run(N_);
  if (factorise_status) {
    return kDecomposerStatusErrorFactorize;
  }

  // save extreme values of the factorisation
  minD_ = factorise.minD_;
  maxD_ = factorise.maxD_;
  minoffD_ = factorise.minoffD_;
  maxoffD_ = factorise.maxoffD_;

  this->valid_ = true;
  use_as_ = false;
  return kDecomposerStatusOk;
}

int FactorHiGHSSolver::solveNE(const HighsSparseMatrix& A,
                               const std::vector<double>& theta,
                               const std::vector<double>& rhs,
                               std::vector<double>& lhs) {
  // only execute the solve if factorization is valid
  assert(this->valid_);

  // initialize lhs with rhs
  lhs = rhs;

  N_.solve(lhs);

  return kDecomposerStatusOk;
}

int FactorHiGHSSolver::solveAS(const HighsSparseMatrix& A,
                               const std::vector<double>& theta,
                               const std::vector<double>& rhs_x,
                               const std::vector<double>& rhs_y,
                               std::vector<double>& lhs_x,
                               std::vector<double>& lhs_y) {
  // only execute the solve if factorization is valid
  assert(this->valid_);

  // create single rhs
  std::vector<double> rhs;
  rhs.insert(rhs.end(), rhs_x.begin(), rhs_x.end());
  rhs.insert(rhs.end(), rhs_y.begin(), rhs_y.end());

  N_.solve(rhs);

  // split lhs
  lhs_x = std::vector<double>(rhs.begin(), rhs.begin() + A.num_col_);
  lhs_y = std::vector<double>(rhs.begin() + A.num_col_, rhs.end());

  return kDecomposerStatusOk;
}

void FactorHiGHSSolver::solveForRefineNE(const HighsSparseMatrix& A,
                                         const std::vector<double>& scaling,
                                         std::vector<double>& res_x,
                                         std::vector<double>& res_y) {
  // size of the matrix
  int n = A.num_col_;
  int m = A.num_row_;

  std::vector<double> temp(n);

  // build rhs = res_y + A Theta res_x
  for (int i = 0; i < n; ++i) {
    temp[i] = res_x[i] / scaling[i];
  }
  A.alphaProductPlusY(1.0, temp, res_y);

  // solve to find cor_y
  N_.solve(res_y);

  // compute cor_x = Theta (A^T cor_y - res_x) into res_x
  A.alphaProductPlusY(-1.0, res_y, res_x, true);
  for (int i = 0; i < n; ++i) {
    res_x[i] /= -scaling[i];
  }
}

void FactorHiGHSSolver::refine(const HighsSparseMatrix& A,
                               const std::vector<double>& scaling,
                               const std::vector<double>& rhs_x,
                               const std::vector<double>& rhs_y,
                               std::vector<double>& delta_x,
                               std::vector<double>& delta_y) {
  // Apply iterative refinement to the augmented system (even if solution was
  // computed with normal equations).

  assert(this->valid_);

  // size of the matrix
  int n = A.num_col_;
  int m = A.num_row_;

  // initialize residual/correction
  std::vector<double> res_x(rhs_x);
  std::vector<double> res_y(rhs_y);

  // residual norm
  double norm_res{};
  double old_norm_res{};
  const double norm_rhs = infNorm(rhs_x, rhs_y);

  // extract dynamic regularizaiton and un-permute it
  std::vector<double> dynamic_reg = N_.dynamic_reg_;
  permuteVector(dynamic_reg, S_.iperm());

  // dynamic_reg corresponds to either primal-dual (for AS), or only dual (for
  // NE). The static primal is already included in theta, so only the static
  // dual needs to be added.

  // compute residual x
  // res_x = rhs_x + (theta^-1 + Rp) * Dx - A^T * Dy
  for (int i = 0; i < n; ++i) {
    // primal dynamic reg is zero for NE. For AS, it is extracted from
    // dynamic_reg.
    double primal_dynamic{};
    if (use_as_) primal_dynamic = dynamic_reg[i];

    res_x[i] += delta_x[i] * (scaling[i] + primal_dynamic);
  }
  A.alphaProductPlusY(-1.0, delta_y, res_x, true);

  //printf("Res x %e\n", infNorm(res_x));

  // compute residual y
  // res_y = rhs_y - A * Dx - Rd * Dy
  A.alphaProductPlusY(-1.0, delta_x, res_y);
  for (int i = 0; i < m; ++i) {
    // dynamir_reg stored only Rd for NE, stores both Rd, Rp for AS
    int offset = use_as_ ? n : 0;
    double dual_reg = kDualStaticRegularization + dynamic_reg[i + offset];
    // if (dynamic_reg[i + offset] > 0.0) printf("%e\n", dual_reg);
    // dual_reg = 0.0;
    res_y[i] -= dual_reg * delta_y[i];
  }

  //printf("Res y %e\n", infNorm(res_y));

  norm_res = infNorm(res_x, res_y);

  for (int iter = 0; iter < 5; ++iter) {
    // stop refinement if residual is small
    if (norm_res < 1e-8) {
      break;
    }

    //printf("%e  --> ", norm_res);

    // compute correction
    // cor = M^-1 res
    if (use_as_) {
      // create single residual
      std::vector<double> res;
      res.insert(res.end(), res_x.begin(), res_x.end());
      res.insert(res.end(), res_y.begin(), res_y.end());

      N_.solve(res);

      // split correction into components
      res_x = std::vector<double>(res.begin(), res.begin() + n);
      res_y = std::vector<double>(res.begin() + n, res.end());
    } else {
      solveForRefineNE(A, scaling, res_x, res_y);
    }

    // add correction to direction
    std::vector<double> temp_x(delta_x);
    std::vector<double> temp_y(delta_y);
    vectorAdd(temp_x, res_x);
    vectorAdd(temp_y, res_y);

    // compute new residual
    res_x = rhs_x;
    for (int i = 0; i < n; ++i) {
      double primal_dynamic{};
      if (use_as_) primal_dynamic = dynamic_reg[i];

      res_x[i] += temp_x[i] * (scaling[i] + primal_dynamic);
    }
    A.alphaProductPlusY(-1.0, temp_y, res_x, true);

    res_y = rhs_y;
    A.alphaProductPlusY(-1.0, temp_x, res_y);
    for (int i = 0; i < m; ++i) {
      int offset = use_as_ ? n : 0;
      double dual_reg = kDualStaticRegularization + dynamic_reg[i + offset];
      res_y[i] -= dual_reg * temp_y[i];
    }

    old_norm_res = norm_res;
    norm_res = infNorm(res_x, res_y);

    // reject if residual became larger
    if (norm_res < old_norm_res) {
      delta_x = temp_x;
      delta_y = temp_y;
    } else {
      //printf(" %e xxx \n", norm_res);
      worst_res_ = std::max(worst_res_, old_norm_res);
      return;
    }
  }

  norm_res = infNorm(res_x, res_y);
  //printf("%e\n", norm_res);

  worst_res_ = std::max(worst_res_, norm_res);
}

void FactorHiGHSSolver::finalise() { S_.printTimes(); }

void FactorHiGHSSolver::extractData(LS_data& data) {
  data.minD = minD_;
  data.maxD = maxD_;
  data.minL = minoffD_;
  data.maxL = maxoffD_;

  // compute number of regularized pivots and maximum regularization
  data.num_reg = 0;
  data.max_reg = 0.0;
  for (double d : N_.dynamic_reg_) {
    if (d != 0.0) {
      ++data.num_reg;
      data.max_reg = std::max(data.max_reg, std::abs(d));
    }
  }

  data.worst_res = worst_res_;
}