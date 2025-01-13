#include "FactorHiGHSSolver.h"

#include "../FactorHiGHS/KrylovMethods.h"

FactorHiGHSSolver::FactorHiGHSSolver(const Options& options)
    : S_((FormatType)options.format), N_(S_) {}

void FactorHiGHSSolver::clear() {
  valid_ = false;
  DataCollector::get()->append();
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
  Analyse analyse(S_, rowsLower, ptrLower, negative_pivots);
  if (int status = analyse.run()) return status;
  DataCollector::get()->printSymbolic(1);

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
  Factorise factorise(S_, rowsLower, ptrLower, valLower);
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
  Factorise factorise(S_, AAt.index_, AAt.start_, AAt.value_);
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
                                         std::vector<double>& res_y) const {
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

void residual_x(const HighsSparseMatrix& A, const std::vector<double>& scaling,
                const std::vector<double>& rhs_x,
                const std::vector<double>& delta_x,
                const std::vector<double>& delta_y,
                const std::vector<double>& total_reg, bool use_as, bool use_reg,
                std::vector<double>& res_x) {
  // compute residual x

  int n = A.num_col_;

  // res_x = rhs_x + (theta^-1 + Rp) * Dx - A^T * Dy
  res_x = rhs_x;
  for (int i = 0; i < n; ++i) {
    // primal reg is only static for NE. For AS, it is extracted from
    // total_reg, already including static and dynamic.
    double primal_reg{};
    if (use_as)
      primal_reg = total_reg[i];
    else
      primal_reg = kPrimalStaticRegularization;

    if (use_reg)
      res_x[i] += delta_x[i] * (scaling[i] + primal_reg);
    else
      res_x[i] += delta_x[i] * scaling[i];
  }
  A.alphaProductPlusY(-1.0, delta_y, res_x, true);
}

void residual_y(const HighsSparseMatrix& A, const std::vector<double>& rhs_y,
                const std::vector<double>& delta_x,
                const std::vector<double>& delta_y,
                const std::vector<double>& total_reg, bool use_as, bool use_reg,
                std::vector<double>& res_y) {
  // compute residual y

  int n = A.num_col_;
  int m = A.num_row_;

  // res_y = rhs_y - A * Dx - Rd * Dy
  res_y = rhs_y;
  A.alphaProductPlusY(-1.0, delta_x, res_y);
  if (use_reg) {
    for (int i = 0; i < m; ++i) {
      // total_reg_ stored only Rd for NE, stores both Rd, Rp for AS
      int offset = use_as ? n : 0;
      double dual_reg = total_reg[i + offset];

      res_y[i] -= dual_reg * delta_y[i];
    }
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
  std::vector<double> total_reg = N_.total_reg_;
  permuteVector(total_reg, S_.iperm());

  // total_reg corresponds to primal-dual (for AS), or only dual (for NE).

  residual_x(A, scaling, rhs_x, delta_x, delta_y, total_reg, use_as_, 1, res_x);
  residual_y(A, rhs_y, delta_x, delta_y, total_reg, use_as_, 1, res_y);
  norm_res = infNorm(res_x, res_y);

#ifdef PRINT_ITER_REF
  printf("Res x %e\n", infNorm(res_x));
  printf("Res y %e\n", infNorm(res_y));
#endif

  for (int iter = 0; iter < kMaxRefinementIter; ++iter) {
    // stop refinement if residual is small
    if (norm_res < kRefinementTolerance) {
      break;
    }

#ifdef PRINT_ITER_REF
    printf("%e  --> ", norm_res);
#endif

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
    residual_x(A, scaling, rhs_x, temp_x, temp_y, total_reg, use_as_, 1, res_x);
    residual_y(A, rhs_y, temp_x, temp_y, total_reg, use_as_, 1, res_y);

    old_norm_res = norm_res;
    norm_res = infNorm(res_x, res_y);

    // reject if residual became larger
    if (norm_res < old_norm_res) {
      delta_x = temp_x;
      delta_y = temp_y;
    } else {
#ifdef PRINT_ITER_REF
      printf(" xxx ");
#endif
      norm_res = old_norm_res;
      break;
    }
  }

#ifdef PRINT_ITER_REF
  printf("%e\n", norm_res);
#endif

  /*{
    IpmMatrix M;
    M.reset(A, scaling, use_as_);
    IpmFactor P(N_);

    if (use_as_) {
      std::vector<double> rhs;
      rhs.insert(rhs.end(), rhs_x.begin(), rhs_x.end());
      rhs.insert(rhs.end(), rhs_y.begin(), rhs_y.end());

      std::vector<double> delta;
      delta.insert(delta.end(), delta_x.begin(), delta_x.end());
      delta.insert(delta.end(), delta_y.begin(), delta_y.end());

      int iter = Gmres(&M, &P, rhs, delta, 1e-6, 100);
      printf("gmres iter %d\n", iter);

      delta_x = std::vector<double>(delta.begin(), delta.begin() + n);
      delta_y = std::vector<double>(delta.begin() + n, delta.end());
    } else {
      std::vector<double> rhs = rhs_y;
      std::vector<double> temp(n);
      for (int i = 0; i < n; ++i) temp[i] = rhs_x[i] / scaling[i];
      A.alphaProductPlusY(1.0, temp, rhs);
      int iter = Cg(&M, &P, rhs, delta_y, 1e-6, 100);
      printf("cg iter %d\n", iter);

      temp = rhs_x;
      A.alphaProductPlusY(-1.0, delta_y, temp, true);
      for (int i = 0; i < n; ++i) delta_x[i] = -temp[i] / scaling[i];
    }
  }*/

  double& p_res = DataCollector::get()->back().perturbed_res;
  double& o_res = DataCollector::get()->back().original_res;

  p_res = std::max(p_res, norm_res);

  // compute residual with respect to original linear system
  residual_x(A, scaling, rhs_x, delta_x, delta_y, total_reg, use_as_, 0, res_x);
  residual_y(A, rhs_y, delta_x, delta_y, total_reg, use_as_, 0, res_y);
  o_res = std::max(o_res, infNorm(res_x, res_y));
}

void FactorHiGHSSolver::finalise() { DataCollector::get()->printTimes(); }