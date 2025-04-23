#include "FactorHiGHSSolver.h"

#include <limits>

namespace highspm {

Int computeLowerAThetaAT(const HighsSparseMatrix& matrix,
                         const std::vector<double>& scaling,
                         HighsSparseMatrix& AAT,
                         const Int max_num_nz = 100000000
                         // Cant exceed kHighsIInf = 2,147,483,647,
                         // otherwise start_ values may overflow. Even
                         // 100,000,000 is probably too large, unless the
                         // matrix is near-full, since fill-in will
                         // overflow pointers
);

FactorHiGHSSolver::FactorHiGHSSolver(const Options& options)
    : S_((FormatType)options.format), N_(S_) {}

void FactorHiGHSSolver::clear() {
  valid_ = false;
  DataCollector::get()->append();
}

Int getNE(const HighsSparseMatrix& A, std::vector<Int>& ptr,
          std::vector<Int>& rows) {
  // Normal equations, full matrix
  std::vector<double> theta;
  HighsSparseMatrix AAt;
  Int status = computeLowerAThetaAT(A, theta, AAt);
  if (status) return kLinearSolverStatusErrorOom;

  rows = std::move(AAt.index_);
  ptr = std::move(AAt.start_);

  return kLinearSolverStatusOk;
}

Int getAS(const HighsSparseMatrix& A, std::vector<Int>& ptr,
          std::vector<Int>& rows) {
  // Augmented system, lower triangular

  Int nA = A.num_col_;
  Int mA = A.num_row_;
  Int nzA = A.numNz();

  ptr.resize(nA + mA + 1);
  rows.resize(nA + nzA + mA);

  Int next = 0;

  for (Int i = 0; i < nA; ++i) {
    // diagonal element
    rows[next] = i;
    ++next;

    // column of A
    for (Int el = A.start_[i]; el < A.start_[i + 1]; ++el) {
      rows[next] = A.index_[el] + nA;
      ++next;
    }

    ptr[i + 1] = next;
  }

  // 2,2 block
  for (Int i = 0; i < mA; ++i) {
    rows[next] = nA + i;
    ++next;
    ptr[nA + i + 1] = ptr[nA + i] + 1;
  }

  return kLinearSolverStatusOk;
}

Int FactorHiGHSSolver::setup(const HighsSparseMatrix& A, Options& options) {
  printf("\n");

  if (Int status = setNla(A, options)) return status;

  setParallel(options);

  DataCollector::get()->printSymbolic(1);
  return kLinearSolverStatusOk;
}

Int FactorHiGHSSolver::factorAS(const HighsSparseMatrix& A,
                                const std::vector<double>& scaling) {
  // only execute factorization if it has not been done yet
  assert(!this->valid_);

  // initialize
  std::vector<Int> ptrLower;
  std::vector<Int> rowsLower;
  std::vector<double> valLower;

  Int nA = A.num_col_;
  Int mA = A.num_row_;
  Int nzA = A.numNz();

  ptrLower.resize(nA + mA + 1);
  rowsLower.resize(nA + nzA + mA);
  valLower.resize(nA + nzA + mA);

  // build lower triangle
  Int next = 0;

  for (Int i = 0; i < nA; ++i) {
    // diagonal element
    rowsLower[next] = i;
    valLower[next++] = -scaling[i];

    // column of A
    for (Int el = A.start_[i]; el < A.start_[i + 1]; ++el) {
      rowsLower[next] = A.index_[el] + nA;
      valLower[next++] = A.value_[el];
    }

    ptrLower[i + 1] = next;
  }

  // 2,2 block
  for (Int i = 0; i < mA; ++i) {
    rowsLower[next] = nA + i;
    valLower[next++] = 0.0;
    ptrLower[nA + i + 1] = ptrLower[nA + i] + 1;
  }

  // factorise matrix
  Factorise factorise(S_, rowsLower, ptrLower, valLower);
  if (factorise.run(N_)) return kLinearSolverStatusErrorFactorise;

  this->valid_ = true;
  use_as_ = true;
  return kLinearSolverStatusOk;
}

Int FactorHiGHSSolver::factorNE(const HighsSparseMatrix& A,
                                const std::vector<double>& scaling) {
  // only execute factorization if it has not been done yet
  assert(!this->valid_);

  // initialize
  std::vector<Int> ptrLower;
  std::vector<Int> rowsLower;
  std::vector<double> valLower;

  Int nA = A.num_col_;
  Int mA = A.num_row_;
  Int nzA = A.numNz();

  // build full matrix
  HighsSparseMatrix AAt;
  Int status = computeLowerAThetaAT(A, scaling, AAt);

  // factorise
  Factorise factorise(S_, AAt.index_, AAt.start_, AAt.value_);
  if (factorise.run(N_)) return kLinearSolverStatusErrorFactorise;

  this->valid_ = true;
  use_as_ = false;
  return kLinearSolverStatusOk;
}

Int FactorHiGHSSolver::solveNE(const std::vector<double>& rhs,
                               std::vector<double>& lhs) {
  // only execute the solve if factorization is valid
  assert(this->valid_);

  // initialize lhs with rhs
  lhs = rhs;

  N_.solve(lhs);

  return kLinearSolverStatusOk;
}

Int FactorHiGHSSolver::solveAS(const std::vector<double>& rhs_x,
                               const std::vector<double>& rhs_y,
                               std::vector<double>& lhs_x,
                               std::vector<double>& lhs_y) {
  // only execute the solve if factorization is valid
  assert(this->valid_);

  Int n = rhs_x.size();

  // create single rhs
  std::vector<double> rhs;
  rhs.insert(rhs.end(), rhs_x.begin(), rhs_x.end());
  rhs.insert(rhs.end(), rhs_y.begin(), rhs_y.end());

  N_.solve(rhs);

  // split lhs
  lhs_x = std::vector<double>(rhs.begin(), rhs.begin() + n);
  lhs_y = std::vector<double>(rhs.begin() + n, rhs.end());

  return kLinearSolverStatusOk;
}

void FactorHiGHSSolver::finalise() { DataCollector::get()->printTimes(); }

Int computeLowerAThetaAT(const HighsSparseMatrix& matrix,
                         const std::vector<double>& scaling,
                         HighsSparseMatrix& AAT, const Int max_num_nz) {
  // Create a row-wise copy of the matrix
  HighsSparseMatrix AT = matrix;
  AT.ensureRowwise();

  Int AAT_dim = matrix.num_row_;
  AAT.num_col_ = AAT_dim;
  AAT.num_row_ = AAT_dim;
  AAT.start_.resize(AAT_dim + 1, 0);

  std::vector<std::tuple<Int, Int, double>> non_zero_values;

  // First pass to calculate the number of non-zero elements in each column
  //
  Int AAT_num_nz = 0;
  std::vector<double> AAT_col_value(AAT_dim, 0);
  std::vector<Int> AAT_col_index(AAT_dim);
  std::vector<bool> AAT_col_in_index(AAT_dim, false);
  for (Int iRow = 0; iRow < AAT_dim; iRow++) {
    // Go along the row of A, and then down the columns corresponding
    // to its nonzeros
    Int num_col_el = 0;
    for (Int iRowEl = AT.start_[iRow]; iRowEl < AT.start_[iRow + 1]; iRowEl++) {
      Int iCol = AT.index_[iRowEl];
      const double theta_value =
          scaling.empty() ? 1.0
                          : 1.0 / (scaling[iCol] + kPrimalStaticRegularization);
      if (!theta_value) continue;
      const double row_value = theta_value * AT.value_[iRowEl];
      for (Int iColEl = matrix.start_[iCol]; iColEl < matrix.start_[iCol + 1];
           iColEl++) {
        Int iRow1 = matrix.index_[iColEl];
        if (iRow1 < iRow) continue;
        double term = row_value * matrix.value_[iColEl];
        if (!AAT_col_in_index[iRow1]) {
          // This entry is not yet in the list of possible nonzeros
          AAT_col_in_index[iRow1] = true;
          AAT_col_index[num_col_el++] = iRow1;
          AAT_col_value[iRow1] = term;
        } else {
          // This entry is in the list of possible nonzeros
          AAT_col_value[iRow1] += term;
        }
      }
    }
    for (Int iEl = 0; iEl < num_col_el; iEl++) {
      Int iCol = AAT_col_index[iEl];
      assert(iCol >= iRow);
      const double value = AAT_col_value[iCol];

      non_zero_values.emplace_back(iRow, iCol, value);
      const Int num_new_nz = 1;
      if (AAT_num_nz + num_new_nz >= max_num_nz)
        return kLinearSolverStatusErrorOom;
      AAT.start_[iRow + 1]++;
      AAT_num_nz += num_new_nz;
      AAT_col_in_index[iCol] = false;
    }
  }

  // Prefix sum to get the correct column pointers
  for (Int i = 0; i < AAT_dim; ++i) AAT.start_[i + 1] += AAT.start_[i];

  AAT.index_.resize(AAT.start_.back());
  AAT.value_.resize(AAT.start_.back());
  AAT.p_end_ = AAT.start_;
  AAT.p_end_.back() = AAT.index_.size();

  std::vector<Int> current_positions = AAT.start_;

  // Second pass to actually fill in the indices and values
  for (const auto& val : non_zero_values) {
    Int i = std::get<0>(val);
    Int j = std::get<1>(val);
    double dot = std::get<2>(val);

    // i>=j, so to get lower triangle, i is the col, j is row
    AAT.index_[current_positions[i]] = j;
    AAT.value_[current_positions[i]] = dot;
    current_positions[i]++;
    AAT.p_end_[i] = current_positions[i];
  }
  AAT.p_end_.clear();
  return kLinearSolverStatusOk;
}

double FactorHiGHSSolver::flops() const { return S_.flops(); }
double FactorHiGHSSolver::spops() const { return S_.spops(); }
double FactorHiGHSSolver::nz() const { return (double)S_.nz(); }

Int FactorHiGHSSolver::choose(const HighsSparseMatrix& A, Options& options) {
  // Choose whether to use augmented system or normal equations.

  assert(options.nla == kOptionNlaChoose);

  Symbolic symb_NE((FormatType)options.format);
  Symbolic symb_AS((FormatType)options.format);
  bool failure_NE = false;
  bool failure_AS = false;

  // Perform analyse phase of normal equations
  {
    std::vector<Int> ptrLower, rowsLower;
    Int NE_status = getNE(A, ptrLower, rowsLower);
    if (NE_status)
      failure_NE = true;
    else {
      Analyse analyse_NE(symb_NE, rowsLower, ptrLower, 0);
      NE_status = analyse_NE.run();

      // NE may fail because of a fail in analyse, or because there are too many
      // nonzeros in the factor, for the given integer type.
      if (NE_status || symb_NE.nz() >= std::numeric_limits<Int>::max())
        failure_NE = true;

      // save data collected for NE and clear record for AS
      DataCollector::get()->saveAndClear();
    }
  }

  // Perform analyse phase of augmented system
  {
    std::vector<Int> ptrLower, rowsLower;
    getAS(A, ptrLower, rowsLower);
    Analyse analyse_AS(symb_AS, rowsLower, ptrLower, A.num_col_);
    Int AS_status = analyse_AS.run();

    // AS may fail because of a fail in analyse, or because there are too many
    // nonzeros in the factor, for the given integer type.
    if (AS_status || symb_AS.nz() >= std::numeric_limits<Int>::max())
      failure_AS = true;
  }

  Int status = kLinearSolverStatusOk;

  // Decision may be forces by failures
  if (failure_NE && !failure_AS) {
    options.nla = kOptionNlaAugmented;
    printf("Using augmented system because normal equations failed\n");
  } else if (failure_AS && !failure_NE) {
    options.nla = kOptionNlaNormEq;
    printf("Using normal equations because augmented system failed\n");
  } else if (failure_AS && failure_NE) {
    status = kLinearSolverStatusErrorAnalyse;
    printf("Failure: both approaches failed analyse phase\n");
  } else {
    // Coefficients for heuristic
    const double kSpopsWeight = 30.0;
    const double kThreshOps = 10.0;
    const double kThreshSn = 1.5;

    // Total number of operations, given by dense flops and sparse indexing
    // operations, weighted with an empirical factor
    double ops_NE = symb_NE.flops() + symb_NE.spops() * kSpopsWeight;
    double ops_AS = symb_AS.flops() + symb_AS.spops() + kSpopsWeight;

    // Average size of supernodes
    double sn_size_NE = (double)symb_NE.size() / symb_NE.sn();
    double sn_size_AS = (double)symb_AS.size() / symb_AS.sn();

    double ratio_ops = ops_NE / ops_AS;
    double ratio_sn = sn_size_AS / sn_size_NE;

    bool NE_much_more_expensive = ratio_ops > kThreshOps;
    bool AS_not_too_expensive = ratio_ops > 1.0 / kThreshOps;
    bool sn_AS_larger_than_NE = ratio_sn > kThreshSn;

    if (NE_much_more_expensive ||
        (sn_AS_larger_than_NE && AS_not_too_expensive)) {
      options.nla = kOptionNlaAugmented;
      printf("Using augmented system because it is preferrable\n");
    } else {
      options.nla = kOptionNlaNormEq;
      printf("Using normal equations because it is preferrable\n");
    }
  }

  if (status != kLinearSolverStatusErrorAnalyse) {
    // DataCollector now contains data of AS in main storage and data of NE in
    // saved storage.

    if (options.nla == kOptionNlaAugmented) {
      S_ = std::move(symb_AS);
      DataCollector::get()->clearSaved();
    } else {
      S_ = std::move(symb_NE);
      DataCollector::get()->loadSaved();
    }
  }

  return status;
}

Int FactorHiGHSSolver::setNla(const HighsSparseMatrix& A, Options& options) {
  std::vector<Int> ptrLower, rowsLower;

  // Build the matrix
  switch (options.nla) {
    case kOptionNlaAugmented: {
      getAS(A, ptrLower, rowsLower);
      Analyse analyse(S_, rowsLower, ptrLower, A.num_col_);
      if (analyse.run()) return kLinearSolverStatusErrorAnalyse;
      printf("Using augmented system as requested\n");
      break;
    }

    case kOptionNlaNormEq: {
      Int NE_status = getNE(A, ptrLower, rowsLower);
      if (NE_status) {
        printf("Failure: AAt is too large\n");
        return kLinearSolverStatusErrorOom;
      }
      Analyse analyse(S_, rowsLower, ptrLower, 0);
      if (analyse.run()) return kLinearSolverStatusErrorAnalyse;
      printf("Using normal equations as requested\n");
      break;
    }

    case kOptionNlaChoose: {
      if (choose(A, options)) return kLinearSolverStatusErrorAnalyse;
      break;
    }
  }

  return kLinearSolverStatusOk;
}

void FactorHiGHSSolver::setParallel(const Options& options) {
  // Set parallel options
  bool parallel_tree = false;
  bool parallel_node = false;

  switch (options.parallel) {
    case kOptionParallelOff:
      printf("Using no parallelism as requested\n");
      break;
    case kOptionParallelOn:
      parallel_tree = true;
      parallel_node = true;
      printf("Using full parallelism as requested\n");
      break;
    case kOptionParallelChoose: {
      parallel_node = true;

      // parallel_tree is active if there is enough parallelism to exploit and
      // if the computational effort is large enough, otherwise the overhead of
      // the scheduler is too much.
      double tree_speedup = S_.flops() / S_.critops();
      if (tree_speedup > kMinTreeSpeedup && S_.flops() > kMinParallelOps) {
        parallel_tree = true;
        printf("Using full parallelism because it is preferrable\n");
      } else {
        printf("Using only node parallelism because it is preferrable\n");
      }
      break;
    }
    case kOptionParallelTreeOnly:
      parallel_tree = true;
      printf("Using only tree parallelism as requested\n");
      break;
    case kOptionParallelNodeOnly:
      parallel_node = true;
      printf("Using only node parallelism as requested\n");
      break;
  }

  S_.setParallel(parallel_tree, parallel_node);
}

}  // namespace highspm