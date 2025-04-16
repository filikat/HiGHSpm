#include "HybridPackedFormatHandler.h"

#include "auxiliary/Auxiliary.h"
#include "CallAndTimeBlas.h"
#include "DataCollector.h"
#include "DenseFact.h"

HybridPackedFormatHandler::HybridPackedFormatHandler(const Symbolic& S, Int sn)
    : FormatHandler(S, sn) {
  // initialize frontal and clique
  initFrontal();
  initClique();
}

void HybridPackedFormatHandler::initFrontal() {
  const Int n_blocks = (sn_size_ - 1) / nb_ + 1;
  diag_start_.resize(n_blocks);
  Int frontal_size = getDiagStart(ldf_, sn_size_, nb_, n_blocks, diag_start_);
  frontal_.resize(frontal_size + extra_space);
  // NB: the plus 10 is not needed, but it avoids weird problems later on.
}

void HybridPackedFormatHandler::initClique() {
  clique_.resize(S_->cliqueSize(sn_));
}

void HybridPackedFormatHandler::assembleFrontal(Int i, Int j, double val) {
  Int block = j / nb_;
  Int ldb = ldf_ - block * nb_;
  Int ii = i - block * nb_;
  Int jj = j - block * nb_;
  frontal_[diag_start_[block] + ii + ldb * jj] = val;
}

void HybridPackedFormatHandler::assembleFrontalMultiple(
    Int num, const std::vector<double>& child, Int nc, Int child_sn, Int row,
    Int col, Int i, Int j) {
  const Int jblock = col / nb_;
  row -= jblock * nb_;
  col -= jblock * nb_;
  const Int start_block = S_->cliqueBlockStart(child_sn, jblock);
  const Int ld = nc - nb_ * jblock;

  Int block = j / nb_;
  Int ldb = ldf_ - block * nb_;
  Int ii = i - block * nb_;
  Int jj = j - block * nb_;

  callAndTime_daxpy(num, 1.0, &child[start_block + row + ld * col], 1,
                    &frontal_[diag_start_[block] + ii + ldb * jj], 1);
}

Int HybridPackedFormatHandler::denseFactorise(double reg_thresh) {
  Int status;

  status = denseFactFP2FH(frontal_.data(), ldf_, sn_size_, nb_);
  if (status) return status;

  // find the position within pivot_sign corresponding to this supernode
  Int sn_start = S_->snStart(sn_);
  const Int* pivot_sign = &S_->pivotSign().data()[sn_start];

  status = denseFactFH(
      'P', ldf_, sn_size_, nb_, frontal_.data(), clique_.data(), pivot_sign,
      reg_thresh, local_reg_.data(), swaps_.data(), pivot_2x2_.data(), sn_);

  return status;
}

void HybridPackedFormatHandler::assembleClique(const std::vector<double>& child,
                                               Int nc, Int child_sn) {
  //   go through the columns of the contribution of the child
  for (Int col = 0; col < nc; ++col) {
    // relative index of column in the frontal matrix
    Int j = S_->relindClique(child_sn, col);

    if (j >= sn_size_) {
      // assemble into clique

      // adjust relative index to access clique
      j -= sn_size_;

      // go through the rows of the contribution of the child
      Int row = col;
      while (row < nc) {
        // relative index of the entry in the matrix clique
        const Int i = S_->relindClique(child_sn, row) - sn_size_;

        // how many entries to sum
        const Int consecutive = S_->consecutiveSums(child_sn, row);

        // use daxpy_ for summing consecutive entries

        const Int jblock_c = col / nb_;
        const Int jb_c = std::min(nb_, nc - nb_ * jblock_c);
        const Int row_c = row - jblock_c * nb_;
        const Int col_c = col - jblock_c * nb_;
        const Int start_block_c = S_->cliqueBlockStart(child_sn, jblock_c);
        const Int ld_c = nc - nb_ * jblock_c;

        const Int jblock = j / nb_;
        const Int jb = std::min(nb_, ldc_ - nb_ * jblock);
        const Int ii = i - jblock * nb_;
        const Int jj = j - jblock * nb_;
        const Int start_block = S_->cliqueBlockStart(sn_, jblock);
        const Int ld = ldc_ - nb_ * jblock;

        callAndTime_daxpy(consecutive, 1.0,
                          &child[start_block_c + row_c + ld_c * col_c], 1,
                          &clique_[start_block + ii + ld * jj], 1);

        row += consecutive;
      }
    }
  }
}

void HybridPackedFormatHandler::extremeEntries() {
  double minD = std::numeric_limits<double>::max();
  double maxD = 0.0;
  double minoffD = std::numeric_limits<double>::max();
  double maxoffD = 0.0;

  // number of blocks of columns
  const Int n_blocks = (sn_size_ - 1) / nb_ + 1;

  // index to access frontal
  Int index{};

  // go through blocks of columns for this supernode
  for (Int j = 0; j < n_blocks; ++j) {
    // number of columns in the block
    const Int jb = std::min(nb_, sn_size_ - nb_ * j);

    for (Int k = 0; k < jb; ++k) {
      // off diagonal entries
      for (Int i = 0; i < k; ++i) {
        if (frontal_[index] != 0.0) {
          minoffD = std::min(minoffD, std::abs(frontal_[index]));
          maxoffD = std::max(maxoffD, std::abs(frontal_[index]));
        }
        index++;
      }

      // diagonal entry
      minD = std::min(minD, std::abs(1.0 / frontal_[index]));
      maxD = std::max(maxD, std::abs(1.0 / frontal_[index]));

      index += jb - k;
    }

    const Int entries_left = (ldf_ - nb_ * j - jb) * jb;

    for (Int i = 0; i < entries_left; ++i) {
      if (frontal_[index] != 0.0) {
        minoffD = std::min(minoffD, std::abs(frontal_[index]));
        maxoffD = std::max(maxoffD, std::abs(frontal_[index]));
      }
      index++;
    }
  }

  DataCollector::get()->setExtremeEntries(minD, maxD, minoffD, maxoffD);
}