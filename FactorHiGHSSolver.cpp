#include "FactorHiGHSSolver.h"

#include "Regularization.h"

FactorHiGHSSolver::FactorHiGHSSolver(const Options& options)
    : S_((FormatType)options.format), N_(S_, data_) {}

void FactorHiGHSSolver::clear() {
  valid_ = false;
  data_.resetExtremeEntries();
}

int FactorHiGHSSolver::setup(const HighsSparseMatrix& A,
                             const Options& options) {
  std::vector<int> ptrLower;
  std::vector<int> rowsLower;

  int nA = A.num_col_;
  int mA = A.num_row_;
  int nzA = A.numNz();

  int negative_pivots{};

  int nla_type = options.nla;

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
  Analyse analyse(S_, data_, rowsLower, ptrLower, negative_pivots);
  if (int status = analyse.run()) return status;
  data_.printSymbolic();

  return kRetOk;
}

int FactorHiGHSSolver::factorAS(const HighsSparseMatrix& A,
                                const std::vector<double>& scaling) {
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
    valLower[next++] = -scaling[i];

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
    valLower[next++] = 0.0;
    ptrLower[nA + i + 1] = ptrLower[nA + i] + 1;
  }

  // factorise matrix
  Factorise factorise(S_, data_, rowsLower, ptrLower, valLower);
  if (factorise.run(N_)) return kDecomposerStatusErrorFactorise;

  this->valid_ = true;
  use_as_ = true;
  return kDecomposerStatusOk;
}

int FactorHiGHSSolver::factorNE(const HighsSparseMatrix& A,
                                const std::vector<double>& scaling) {
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
  int status = computeLowerAThetaAT(A, scaling, AAt);

  // factorise
  Factorise factorise(S_, data_, AAt.index_, AAt.start_, AAt.value_);
  if (factorise.run(N_)) return kDecomposerStatusErrorFactorise;

  this->valid_ = true;
  use_as_ = false;
  return kDecomposerStatusOk;
}

int FactorHiGHSSolver::solveNE(const std::vector<double>& rhs,
                               std::vector<double>& lhs) {
  // only execute the solve if factorization is valid
  assert(this->valid_);

  // initialize lhs with rhs
  lhs = rhs;

  N_.solve(lhs);

  return kDecomposerStatusOk;
}

int FactorHiGHSSolver::solveAS(const std::vector<double>& rhs_x,
                               const std::vector<double>& rhs_y,
                               std::vector<double>& lhs_x,
                               std::vector<double>& lhs_y) {
  // only execute the solve if factorization is valid
  assert(this->valid_);

  int n = rhs_x.size();

  // create single rhs
  std::vector<double> rhs;
  rhs.insert(rhs.end(), rhs_x.begin(), rhs_x.end());
  rhs.insert(rhs.end(), rhs_y.begin(), rhs_y.end());

  N_.solve(rhs);

  // split lhs
  lhs_x = std::vector<double>(rhs.begin(), rhs.begin() + n);
  lhs_y = std::vector<double>(rhs.begin() + n, rhs.end());

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

  // build rhs = res_y + A (Theta^-1+Rp)^-1 res_x
  for (int i = 0; i < n; ++i) {
    temp[i] = res_x[i] / (scaling[i] + kPrimalStaticRegularization);
  }
  A.alphaProductPlusY(1.0, temp, res_y);

  // solve to find cor_y
  N_.solve(res_y);

  // compute cor_x = (Theta^-1+Rp)^-1 (A^T cor_y - res_x) into res_x
  A.alphaProductPlusY(-1.0, res_y, res_x, true);
  for (int i = 0; i < n; ++i) {
    res_x[i] /= -(scaling[i] + kPrimalStaticRegularization);
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
  std::vector<double> res_x;
  std::vector<double> res_y;

  // residual norm
  double norm_res{};
  double old_norm_res{};
  const double norm_rhs = infNorm(rhs_x, rhs_y);

  // extract dynamic regularizaiton and un-permute it
  std::vector<double> total_reg_ = N_.total_reg_;
  permuteVector(total_reg_, S_.iperm());

  // dynamic_reg corresponds to either primal-dual (for AS), or only dual (for
  // NE).

  // compute residual x
  // res_x = rhs_x + (theta^-1 + Rp) * Dx - A^T * Dy
  res_x = rhs_x;
  for (int i = 0; i < n; ++i) {
    // primal reg is only static for NE. For AS, it is extracted from
    // total_reg, already including static and dynamic.
    double primal_reg{};
    if (use_as_)
      primal_reg = total_reg_[i];
    else
      primal_reg = kPrimalStaticRegularization;

    res_x[i] += delta_x[i] * (scaling[i] + primal_reg);
  }
  A.alphaProductPlusY(-1.0, delta_y, res_x, true);

  // printf("Res x %e\n", infNorm(res_x));

  // compute residual y
  // res_y = rhs_y - A * Dx - Rd * Dy
  res_y = rhs_y;
  A.alphaProductPlusY(-1.0, delta_x, res_y);
  for (int i = 0; i < m; ++i) {
    // total_reg_ stored only Rd for NE, stores both Rd, Rp for AS
    int offset = use_as_ ? n : 0;
    double dual_reg = total_reg_[i + offset];

    res_y[i] -= dual_reg * delta_y[i];
  }

  // printf("Res y %e\n", infNorm(res_y));

  norm_res = infNorm(res_x, res_y);

  for (int iter = 0; iter < 5; ++iter) {
    // stop refinement if residual is small
    if (norm_res < 1e-8) {
      break;
    }

    // printf("%e  --> ", norm_res);

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
      double primal_reg{};
      if (use_as_)
        primal_reg = total_reg_[i];
      else
        primal_reg = kPrimalStaticRegularization;

      res_x[i] += temp_x[i] * (scaling[i] + primal_reg);
    }
    A.alphaProductPlusY(-1.0, temp_y, res_x, true);

    res_y = rhs_y;
    A.alphaProductPlusY(-1.0, temp_x, res_y);
    for (int i = 0; i < m; ++i) {
      int offset = use_as_ ? n : 0;
      double dual_reg = total_reg_[i + offset];
      res_y[i] -= dual_reg * temp_y[i];
    }

    old_norm_res = norm_res;
    norm_res = infNorm(res_x, res_y);

    // reject if residual became larger
    if (norm_res < old_norm_res) {
      delta_x = temp_x;
      delta_y = temp_y;
    } else {
      // printf(" %e xxx \n", norm_res);
      data_.setWorstRes(old_norm_res);
      return;
    }
  }

  norm_res = infNorm(res_x, res_y);
  // printf("%e\n", norm_res);

  data_.setWorstRes(norm_res);
}

void FactorHiGHSSolver::finalise() { data_.printTimes(); }