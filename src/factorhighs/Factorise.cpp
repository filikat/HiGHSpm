#include "Factorise.h"

#include <algorithm>
#include <fstream>

#include "auxiliary/Auxiliary.h"
#include "DataCollector.h"
#include "FactorHiGHSSettings.h"
#include "FormatHandler.h"
#include "FullFormatHandler.h"
#include "HybridHybridFormatHandler.h"
#include "HybridPackedFormatHandler.h"
#include "PackedPackedFormatHandler.h"
#include "ReturnValues.h"
#include "SymScaling.h"
#include "parallel/HighsParallel.h"

Factorise::Factorise(const Symbolic& S, const std::vector<int>& rowsA,
                     const std::vector<int>& ptrA,
                     const std::vector<double>& valA)
    : S_{S} {
  // Input the symmetric matrix to be factorised in CSC format and the symbolic
  // factorisation coming from Analyze.
  // Only the lower triangular part of the matrix is used.
  // The Factorise object takes ownership of the matrix; rowsA, ptrA and valA
  // are not valid anymore.

  n_ = ptrA.size() - 1;

  if (n_ != S_.size()) {
    printf(
        "Matrix provided to Factorise has size incompatible with symbolic "
        "object.\n");
    return;
  }

  // take ownership of the matrix
  rowsA_ = std::move(rowsA);
  valA_ = std::move(valA);
  ptrA_ = std::move(ptrA);

  // Permute the matrix.
  // This also removes any entry not in the lower triangle.
  permute(S_.iperm());

  nzA_ = ptrA_.back();

  // Double transpose to sort columns
  std::vector<int> temp_ptr(n_ + 1);
  std::vector<int> temp_rows(nzA_);
  std::vector<double> temp_val(nzA_);
  transpose(ptrA_, rowsA_, valA_, temp_ptr, temp_rows, temp_val);
  transpose(temp_ptr, temp_rows, temp_val, ptrA_, rowsA_, valA_);

  // create linked lists of children in supernodal elimination tree
  childrenLinkedList(S_.snParent(), first_child_, next_child_);

#ifdef PARALLEL_TREE
  // create reverse linked lists of children
  first_child_reverse_ = first_child_;
  next_child_reverse_ = next_child_;
  reverseLinkedList(first_child_reverse_, next_child_reverse_);
#endif

  // compute largest diagonal entry in absolute value
  max_diag_ = 0.0;
  min_diag_ = std::numeric_limits<double>::max();
  for (int col = 0; col < n_; ++col) {
    double val = std::abs(valA_[ptrA_[col]]);
    max_diag_ = std::max(max_diag_, val);
    min_diag_ = std::min(min_diag_, val);
  }

  // compute norm 1 of matrix
  std::vector<double> col_norm1(n_);
  for (int col = 0; col < n_; ++col) {
    for (int el = ptrA_[col]; el < ptrA_[col + 1]; ++el) {
      int row = rowsA_[el];
      double val = valA_[el];
      col_norm1[col] += std::abs(val);
      if (row != col) col_norm1[row] += std::abs(val);
    }
  }
  A_norm1_ = *std::max_element(col_norm1.begin(), col_norm1.end());

  DataCollector::get()->back().M_norm1 = A_norm1_;
  DataCollector::get()->back().M_maxdiag = max_diag_;

  /*print(S_.snParent(), "parent");
  print(S_.snStart(), "start");
  print(rowsA_, "Mrows");
  print(ptrA_, "Mptr");
  print(valA_, "Mval");*/
}

void Factorise::permute(const std::vector<int>& iperm) {
  // Symmetric permutation of the lower triangular matrix A based on inverse
  // permutation iperm.
  // The resulting matrix is lower triangular, regardless of the input matrix.

  std::vector<int> work(n_, 0);

  // go through the columns to count the nonzeros
  for (int j = 0; j < n_; ++j) {
    // get new index of column
    const int col = iperm[j];

    // go through elements of column
    for (int el = ptrA_[j]; el < ptrA_[j + 1]; ++el) {
      const int i = rowsA_[el];

      // ignore potential entries in upper triangular part
      if (i < j) continue;

      // get new index of row
      const int row = iperm[i];

      // since only lower triangular part is used, col is smaller than row
      int actual_col = std::min(row, col);
      ++work[actual_col];
    }
  }

  std::vector<int> new_ptr(n_ + 1);

  // get column pointers by summing the count of nonzeros in each column.
  // copy column pointers into work
  counts2Ptr(new_ptr, work);

  std::vector<int> new_rows(new_ptr.back());
  std::vector<double> new_val(new_ptr.back());

  // go through the columns to assign row indices
  for (int j = 0; j < n_; ++j) {
    // get new index of column
    const int col = iperm[j];

    // go through elements of column
    for (int el = ptrA_[j]; el < ptrA_[j + 1]; ++el) {
      const int i = rowsA_[el];

      // ignore potential entries in upper triangular part
      if (i < j) continue;

      // get new index of row
      const int row = iperm[i];

      // since only lower triangular part is used, col is smaller than row
      const int actual_col = std::min(row, col);
      const int actual_row = std::max(row, col);

      int pos = work[actual_col]++;
      new_rows[pos] = actual_row;
      new_val[pos] = valA_[el];
    }
  }

  ptrA_ = std::move(new_ptr);
  rowsA_ = std::move(new_rows);
  valA_ = std::move(new_val);
}

std::unique_ptr<FormatHandler> getFormatHandler(const Symbolic& S, int sn) {
  std::unique_ptr<FormatHandler> ptr;
  switch (S.formatType()) {
    case FormatType::Full:
      ptr.reset(new FullFormatHandler(S, sn));
      break;
    case FormatType::HybridPacked:
      ptr.reset(new HybridPackedFormatHandler(S, sn));
      break;
    case FormatType::HybridHybrid:
      ptr.reset(new HybridHybridFormatHandler(S, sn));
      break;
    case FormatType::PackedPacked:
      ptr.reset(new PackedPackedFormatHandler(S, sn));
      break;
  }
  return ptr;
}

void Factorise::processSupernode(int sn) {
  // Assemble frontal matrix for supernode sn, perform partial factorisation and
  // store the result.

  if (flag_stop_) return;

#ifdef PARALLEL_TREE
  // thr_per_sn[sn] = highs::parallel::thread_num();

  // spawn children of this supernode in reverse order
  int child_to_spawn = first_child_reverse_[sn];
  while (child_to_spawn != -1) {
    highs::parallel::spawn([=]() { processSupernode(child_to_spawn); });
    child_to_spawn = next_child_reverse_[child_to_spawn];
  }

  // wait for first child to finish, before starting the parent (if there is a
  // first child)
  if (first_child_reverse_[sn] != -1) highs::parallel::sync();
#endif

#ifdef FINE_TIMING
  Clock clock;
#endif
  // ===================================================
  // Supernode information
  // ===================================================
  // first and last+1 column of the supernodes
  const int sn_begin = S_.snStart(sn);
  const int sn_end = S_.snStart(sn + 1);
  const int sn_size = sn_end - sn_begin;

  // initialize the format handler
  // this also allocates space for the frontal matrix and schur complement
  std::unique_ptr<FormatHandler> FH = getFormatHandler(S_, sn);

#ifdef FINE_TIMING
  DataCollector::get()->sumTime(kTimeFactorisePrepare, clock.stop());
#endif

#ifdef FINE_TIMING
  clock.start();
#endif
  // ===================================================
  // Assemble original matrix A into frontal
  // ===================================================
  // j is relative column index in the frontal matrix
  for (int j = 0; j < sn_size; ++j) {
    // column index in the original matrix
    const int col = sn_begin + j;

    // go through the column
    for (int el = ptrA_[col]; el < ptrA_[col + 1]; ++el) {
      // relative row index in the frontal matrix
      const int i = S_.relindCols(el);

      FH->assembleFrontal(i, j, valA_[el]);
    }
  }
#ifdef FINE_TIMING
  DataCollector::get()->sumTime(kTimeFactoriseAssembleOriginal, clock.stop());
#endif

  // ===================================================
  // Assemble frontal matrices of children
  // ===================================================
  int child_sn = first_child_[sn];
  while (child_sn != -1) {
    // Schur contribution of the current child
    std::vector<double>& child_clique = schur_contribution_[child_sn];

#ifdef PARALLEL_TREE
    // sync with spawned child, apart from the first one
    if (child_sn != first_child_[sn]) highs::parallel::sync();

    if (flag_stop_) return;

    if (child_clique.size() == 0) {
      printf("Missing child supernode contribution\n");
      flag_stop_ = true;
      return;
    }
#endif

    // determine size of clique of child
    const int child_begin = S_.snStart(child_sn);
    const int child_end = S_.snStart(child_sn + 1);

    // number of nodes in child sn
    const int child_size = child_end - child_begin;

    // size of clique of child sn
    const int nc = S_.ptr(child_sn + 1) - S_.ptr(child_sn) - child_size;

// ASSEMBLE INTO FRONTAL
#ifdef FINE_TIMING
    clock.start();
#endif
    // go through the columns of the contribution of the child
    for (int col = 0; col < nc; ++col) {
      // relative index of column in the frontal matrix
      int j = S_.relindClique(child_sn, col);

      if (j < sn_size) {
        // assemble into frontal

        // go through the rows of the contribution of the child
        int row = col;
        while (row < nc) {
          // relative index of the entry in the matrix frontal
          const int i = S_.relindClique(child_sn, row);

          // how many entries to sum
          const int consecutive = S_.consecutiveSums(child_sn, row);

          FH->assembleFrontalMultiple(consecutive, child_clique, nc, child_sn,
                                      row, col, i, j);

          row += consecutive;
        }
      }
    }
#ifdef FINE_TIMING
    DataCollector::get()->sumTime(kTimeFactoriseAssembleChildrenFrontal,
                                  clock.stop());
#endif

// ASSEMBLE INTO CLIQUE
#ifdef FINE_TIMING
    clock.start();
#endif
    FH->assembleClique(child_clique, nc, child_sn);
#ifdef FINE_TIMING
    DataCollector::get()->sumTime(kTimeFactoriseAssembleChildrenClique,
                                  clock.stop());
#endif

    // Schur contribution of the child is no longer needed
    // Swap with temporary empty vector to deallocate memory
    std::vector<double> temp_empty;
    schur_contribution_[child_sn].swap(temp_empty);

    // move on to the next child
    child_sn = next_child_[child_sn];
  }

  if (flag_stop_) return;

    // ===================================================
    // Partial factorisation
    // ===================================================
#ifdef FINE_TIMING
  clock.start();
#endif
  // threshold for regularization
  // const double reg_thresh = max_diag_ * kDynamicDiagCoeff;
  const double reg_thresh = A_norm1_ * kDynamicDiagCoeff;

  if (FH->denseFactorise(reg_thresh)) {
    flag_stop_ = true;
    return;
  }
#ifdef FINE_TIMING
  DataCollector::get()->sumTime(kTimeFactoriseDenseFact, clock.stop());
#endif

#ifdef FINE_TIMING
  clock.start();
#endif
  // compute largest elements in factorization
  FH->extremeEntries();

  // terminate the format handler
  FH->terminate(sn_columns_[sn], schur_contribution_[sn], total_reg_,
                swaps_[sn], pivot_2x2_[sn]);
#ifdef FINE_TIMING
  DataCollector::get()->sumTime(kTimeFactoriseTerminate, clock.stop());
#endif
}

bool Factorise::run(Numeric& num) {
#ifdef COARSE_TIMING
  Clock clock;
#endif

  total_reg_.assign(n_, 0.0);

  // allocate space for list of generated elements and columns of L
  schur_contribution_.resize(S_.sn());
  sn_columns_.resize(S_.sn());
  swaps_.resize(S_.sn());
  pivot_2x2_.resize(S_.sn());

#ifdef PARALLEL_TREE
  int spawned_roots{};
  // spawn tasks for root supernodes
  for (int sn = 0; sn < S_.sn(); ++sn) {
    if (S_.snParent(sn) == -1) {
      highs::parallel::spawn([=]() { processSupernode(sn); });
      ++spawned_roots;
    }
  }

  // sync tasks for root supernodes
  for (int root = 0; root < spawned_roots; ++root) {
    highs::parallel::sync();
  }

#else
  // go through each supernode serially
  for (int sn = 0; sn < S_.sn(); ++sn) {
    processSupernode(sn);
  }
#endif

  if (flag_stop_) return true;

  // move factorisation to numerical object
  num.sn_columns_ = std::move(sn_columns_);
  num.total_reg_ = std::move(total_reg_);
  num.swaps_ = std::move(swaps_);
  num.pivot_2x2_ = std::move(pivot_2x2_);
  num.ptrA_ = std::move(ptrA_);
  num.rowsA_ = std::move(rowsA_);
  num.valA_ = std::move(valA_);

#ifdef COARSE_TIMING
  DataCollector::get()->sumTime(kTimeFactorise, clock.stop());
#endif

  return false;
}
