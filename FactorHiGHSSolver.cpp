#include "FactorHiGHSSolver.h"

#include "../FactorHiGHS/KrylovMethods.h"

int computeLowerAThetaAT(const HighsSparseMatrix& matrix,
                         const std::vector<double>& scaling,
                         HighsSparseMatrix& AAT,
                         const int max_num_nz = 100000000
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
      return kLinearSolverStatusErrorOom;
    }

    rowsLower = std::move(AAt.index_);
    ptrLower = std::move(AAt.start_);
  }

  // Perform analyse phase
  Analyse analyse(S_, rowsLower, ptrLower, negative_pivots);
  if (int status = analyse.run()) return kLinearSolverStatusErrorAnalyse;
  DataCollector::get()->printSymbolic(1);

  return kLinearSolverStatusOk;
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
  if (factorise.run(N_)) return kLinearSolverStatusErrorFactorise;

  this->valid_ = true;
  use_as_ = true;
  return kLinearSolverStatusOk;
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
  if (factorise.run(N_)) return kLinearSolverStatusErrorFactorise;

  this->valid_ = true;
  use_as_ = false;
  return kLinearSolverStatusOk;
}

int FactorHiGHSSolver::solveNE(const std::vector<double>& rhs,
                               std::vector<double>& lhs) {
  // only execute the solve if factorization is valid
  assert(this->valid_);

  // initialize lhs with rhs
  lhs = rhs;

  N_.solve(lhs);

  return kLinearSolverStatusOk;
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

  return kLinearSolverStatusOk;
}

void FactorHiGHSSolver::finalise() { DataCollector::get()->printTimes(); }

int computeLowerAThetaAT(const HighsSparseMatrix& matrix,
                         const std::vector<double>& scaling,
                         HighsSparseMatrix& AAT, const int max_num_nz) {
  // Create a row-wise copy of the matrix
  HighsSparseMatrix AT = matrix;
  AT.ensureRowwise();

  int AAT_dim = matrix.num_row_;
  AAT.num_col_ = AAT_dim;
  AAT.num_row_ = AAT_dim;
  AAT.start_.resize(AAT_dim + 1, 0);

  std::vector<std::tuple<int, int, double>> non_zero_values;

  // First pass to calculate the number of non-zero elements in each column
  //
  int AAT_num_nz = 0;
  std::vector<double> AAT_col_value(AAT_dim, 0);
  std::vector<int> AAT_col_index(AAT_dim);
  std::vector<bool> AAT_col_in_index(AAT_dim, false);
  for (int iRow = 0; iRow < AAT_dim; iRow++) {
    // Go along the row of A, and then down the columns corresponding
    // to its nonzeros
    int num_col_el = 0;
    for (int iRowEl = AT.start_[iRow]; iRowEl < AT.start_[iRow + 1]; iRowEl++) {
      int iCol = AT.index_[iRowEl];
      const double theta_value =
          scaling.empty() ? 1.0
                          : 1.0 / (scaling[iCol] + kPrimalStaticRegularization);
      if (!theta_value) continue;
      const double row_value = theta_value * AT.value_[iRowEl];
      for (int iColEl = matrix.start_[iCol]; iColEl < matrix.start_[iCol + 1];
           iColEl++) {
        int iRow1 = matrix.index_[iColEl];
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
    for (int iEl = 0; iEl < num_col_el; iEl++) {
      int iCol = AAT_col_index[iEl];
      assert(iCol >= iRow);
      const double value = AAT_col_value[iCol];

      non_zero_values.emplace_back(iRow, iCol, value);
      const int num_new_nz = 1;
      if (AAT_num_nz + num_new_nz >= max_num_nz)
        return kLinearSolverStatusErrorOom;
      AAT.start_[iRow + 1]++;
      AAT_num_nz += num_new_nz;
      AAT_col_in_index[iCol] = false;
    }
  }

  // Prefix sum to get the correct column pointers
  for (int i = 0; i < AAT_dim; ++i) AAT.start_[i + 1] += AAT.start_[i];

  AAT.index_.resize(AAT.start_.back());
  AAT.value_.resize(AAT.start_.back());
  AAT.p_end_ = AAT.start_;
  AAT.p_end_.back() = AAT.index_.size();

  std::vector<int> current_positions = AAT.start_;

  // Second pass to actually fill in the indices and values
  for (const auto& val : non_zero_values) {
    int i = std::get<0>(val);
    int j = std::get<1>(val);
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
double FactorHiGHSSolver::nz() const { return S_.nz(); }