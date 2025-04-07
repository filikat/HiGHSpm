#include "FullFormatHandler.h"

#include "CallAndTimeBlas.h"
#include "DataCollector.h"
#include "DenseFact.h"

FullFormatHandler::FullFormatHandler(const Symbolic& S, int sn)
    : FormatHandler(S, sn) {
  // initialize frontal and clique
  initFrontal();
  initClique();
}

void FullFormatHandler::initFrontal() {
  // frontal is initialized to zero
  frontal_.resize(ldf_ * sn_size_);
}

void FullFormatHandler::initClique() { return clique_.resize(ldc_ * ldc_); }

void FullFormatHandler::assembleFrontal(int i, int j, double val) {
  frontal_[i + j * ldf_] = val;
}

void FullFormatHandler::assembleFrontalMultiple(
    int num, const std::vector<double>& child, int nc, int child_sn, int row,
    int col, int i, int j) {
  callAndTime_daxpy(num, 1.0, &child[row + nc * col], 1,
                    &frontal_[i + ldf_ * j], 1);
}

int FullFormatHandler::denseFactorise(double reg_thresh) {
  int status;

  // find the position within pivot_sign corresponding to this supernode
  int sn_start = S_->snStart(sn_);
  const int* pivot_sign = &S_->pivotSign().data()[sn_start];

  status =
      denseFactF(ldf_, sn_size_, nb_, frontal_.data(), ldf_, clique_.data(),
                 ldc_, pivot_sign, reg_thresh, local_reg_.data(), sn_);

  return status;
}

void FullFormatHandler::assembleClique(const std::vector<double>& child, int nc,
                                       int child_sn) {
  //   go through the columns of the contribution of the child
  for (int col = 0; col < nc; ++col) {
    // relative index of column in the frontal matrix
    int j = S_->relindClique(child_sn, col);

    if (j >= sn_size_) {
      // assemble into clique

      // adjust relative index to access clique
      j -= sn_size_;

      // go through the rows of the contribution of the child
      int row = col;
      while (row < nc) {
        // relative index of the entry in the matrix clique
        const int i = S_->relindClique(child_sn, row) - sn_size_;

        // how many entries to sum
        const int consecutive = S_->consecutiveSums(child_sn, row);

        // use daxpy_ for summing consecutive entries
        callAndTime_daxpy(consecutive, 1.0, &child[row + nc * col], 1,
                          &clique_[i + ldc_ * j], 1);

        row += consecutive;
      }
    }
  }
}

void FullFormatHandler::extremeEntries() {
  double minD = 1e100;
  double maxD = 0.0;
  double minoffD = 1e100;
  double maxoffD = 0.0;

  for (int col = 0; col < sn_size_; ++col) {
    // diagonal entry
    minD = std::min(minD, std::abs(1.0 / frontal_[col + ldf_ * col]));
    maxD = std::max(maxD, std::abs(1.0 / frontal_[col + ldf_ * col]));

    // off diagonal entries
    for (int row = col + 1; row < ldf_; ++row) {
      if (frontal_[row + ldf_ * col] != 0.0) {
        minoffD = std::min(minoffD, std::abs(frontal_[row + ldf_ * col]));
        maxoffD = std::max(maxoffD, std::abs(frontal_[row + ldf_ * col]));
      }
    }
  }

  DataCollector::get()->setExtremeEntries(minD, maxD, minoffD, maxoffD);
}