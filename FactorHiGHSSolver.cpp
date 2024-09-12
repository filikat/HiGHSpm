#include "FactorHiGHSSolver.h"

#include "MA86Solver.h"

void FactorHiGHSSolver::clear() { valid_ = false; }

void FactorHiGHSSolver::setup(const HighsSparseMatrix& A, int type) {
  std::vector<int> ptrLower;
  std::vector<int> rowsLower;

  int nA = A.num_col_;
  int mA = A.num_row_;
  int nzA = A.numNz();

  FactType fact_type = (type == 0) ? AugSys : NormEq;

  // Build the matrix
  if (fact_type == FactType::AugSys) {
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

  } else {
    // Normal equations, full matrix
    std::vector<double> theta;
    HighsSparseMatrix AAt;
    int status = computeAThetaAT(A, theta, AAt);

    // extract lower triangle
    for (int col = 0; col < mA; ++col) {
      ptrLower.push_back(rowsLower.size());
      for (int el = AAt.start_[col]; el < AAt.start_[col + 1]; el++) {
        int row = AAt.index_[el];
        if (row >= col) {
          rowsLower.push_back(row);
        }
      }
    }
    ptrLower.push_back(rowsLower.size());
  }

  // Perform analyse phase
  Analyse analyse(rowsLower, ptrLower, fact_type);
  analyse.run(S_);
  //S_.print();
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
    valLower[next++] = 0.0;
    ptrLower[nA + i + 1] = ptrLower[nA + i] + 1;
  }

  // factorise matrix
  Factorise factorise(S_, rowsLower, ptrLower, valLower);
  int status = factorise.run(N_);
  if (status) {
    return kDecomposerStatusErrorFactorize;
  }

  this->valid_ = true;
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
  int status = computeAThetaAT(A, theta, AAt);

  // extract lower triangle
  for (int col = 0; col < mA; ++col) {
    ptrLower.push_back(valLower.size());
    for (int el = AAt.start_[col]; el < AAt.start_[col + 1]; el++) {
      int row = AAt.index_[el];
      if (row >= col) {
        valLower.push_back(AAt.value_[el]);
        rowsLower.push_back(row);
      }
    }
  }
  ptrLower.push_back(valLower.size());

  // factorise
  Factorise factorise(S_, rowsLower, ptrLower, valLower);
  int factorise_status = factorise.run(N_);
  if (factorise_status) {
    return kDecomposerStatusErrorFactorize;
  }

  this->valid_ = true;
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
