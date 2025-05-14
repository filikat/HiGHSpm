#include "IpmMatrix.h"

#include "factorhighs/FactorHiGHSSettings.h"

namespace highspm {

IpmMatrix::IpmMatrix(const HighsSparseMatrix& A,
                     const std::vector<double>& scaling, bool use_as)
    : A_{&A},
      scaling_{scaling.data()},
      n_{(Int)scaling.size()},
      use_as_{use_as} {}

void IpmMatrix::applyNE(std::vector<double>& x) const {
  // compute A * Theta * A^T * x

  std::vector<double> temp(n_);

  // temp = A^T * x
  A_->alphaProductPlusY(1.0, x, temp, true);

  // temp = Theta * temp
  for (Int i = 0; i < n_; ++i)
    temp[i] /= (scaling_[i] + kPrimalStaticRegularisation);

  // x = A * temp
  x.assign(x.size(), 0.0);
  A_->alphaProductPlusY(1.0, temp, x);
}

void IpmMatrix::applyNE(std::vector<HighsCDouble>& x) const {
  // compute A * Theta * A^T * x

  std::vector<HighsCDouble> temp(n_);

  // HighsSparseMatrix does not have a function to do products with
  // HighsCDouble, so do it manually.

  // temp = A^T * x
  for (Int iCol = 0; iCol < A_->num_col_; iCol++)
    for (Int iEl = A_->start_[iCol]; iEl < A_->start_[iCol + 1]; iEl++)
      temp[iCol] += (HighsCDouble)A_->value_[iEl] * x[A_->index_[iEl]];

  // temp = Theta * temp
  for (Int i = 0; i < n_; ++i)
    temp[i] /= (HighsCDouble)(scaling_[i] + kPrimalStaticRegularisation);

  // x = A * temp
  x.assign(x.size(), 0.0);
  for (Int iCol = 0; iCol < A_->num_col_; iCol++)
    for (Int iEl = A_->start_[iCol]; iEl < A_->start_[iCol + 1]; iEl++)
      x[A_->index_[iEl]] += (HighsCDouble)A_->value_[iEl] * temp[iCol];
}

void IpmMatrix::applyAS(std::vector<double>& x) const {
  // compute y1 = -Theta^-1 * x1 + A^T * x2
  //         y2 = A * x1

  std::vector<double> x1(x.begin(), x.begin() + n_);
  std::vector<double> x2(x.begin() + n_, x.end());

  std::vector<double> y1(A_->num_col_);
  std::vector<double> y2(A_->num_row_);

  for (Int i = 0; i < n_; ++i) y1[i] = -x1[i] * scaling_[i];
  A_->alphaProductPlusY(1.0, x2, y1, true);

  A_->alphaProductPlusY(1.0, x1, y2);

  std::memcpy(x.data(), y1.data(), y1.size() * sizeof(double));
  std::memcpy(x.data() + n_, y2.data(), y2.size() * sizeof(double));
}

void IpmMatrix::applyAS(std::vector<HighsCDouble>& x) const {
  // compute y1 = -Theta^-1 * x1 + A^T * x2
  //         y2 = A * x1

  std::vector<HighsCDouble> x1(x.begin(), x.begin() + n_);
  std::vector<HighsCDouble> x2(x.begin() + n_, x.end());

  x.assign(x.size(), 0.0);

  HighsCDouble* y1 = x.data();
  HighsCDouble* y2 = x.data() + n_;

  for (Int i = 0; i < n_; ++i) y1[i] = -x1[i] * (HighsCDouble)scaling_[i];

  for (Int iCol = 0; iCol < A_->num_col_; iCol++)
    for (Int iEl = A_->start_[iCol]; iEl < A_->start_[iCol + 1]; iEl++)
      y1[iCol] += (HighsCDouble)A_->value_[iEl] * x2[A_->index_[iEl]];

  for (Int iCol = 0; iCol < A_->num_col_; iCol++)
    for (Int iEl = A_->start_[iCol]; iEl < A_->start_[iCol + 1]; iEl++)
      y2[A_->index_[iEl]] += (HighsCDouble)A_->value_[iEl] * x1[iCol];
}

void IpmMatrix::apply(std::vector<double>& x) const {
  if (use_as_)
    applyAS(x);
  else
    applyNE(x);
}
void IpmMatrix::apply(std::vector<HighsCDouble>& x) const {
  if (use_as_)
    applyAS(x);
  else
    applyNE(x);
}

}  // namespace highspm