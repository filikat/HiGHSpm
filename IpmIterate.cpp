#include "IpmIterate.h"

#include "../FactorHiGHS/DataCollector.h"

IpmIterate::IpmIterate(const IpmModel& model) : model_{model} {
  x_.resize(model_.n_);
  xl_.resize(model_.n_);
  xu_.resize(model_.n_);
  y_.resize(model_.m_);
  zl_.resize(model_.n_);
  zu_.resize(model_.n_);

  clearRes();
}

bool IpmIterate::isNan() const {
  if (isNanVector(x_) || isNanVector(xl_) || isNanVector(xu_) ||
      isNanVector(y_) || isNanVector(zl_) || isNanVector(zu_) ||
      isNanVector(res1_) || isNanVector(res2_) || isNanVector(res3_) ||
      isNanVector(res4_) || isNanVector(res5_) || isNanVector(res6_))
    return true;
  return false;
}

bool IpmIterate::isInf() const {
  if (isInfVector(x_) || isInfVector(xl_) || isInfVector(xu_) ||
      isInfVector(y_) || isInfVector(zl_) || isInfVector(zu_) ||
      isInfVector(res1_) || isInfVector(res2_) || isInfVector(res3_) ||
      isInfVector(res4_) || isInfVector(res5_) || isInfVector(res6_))
    return true;
  return false;
}

void IpmIterate::mu() {
  mu_ = 0.0;
  int number_finite_bounds{};
  for (int i = 0; i < model_.n_; ++i) {
    if (model_.hasLb(i)) {
      mu_ += xl_[i] * zl_[i];
      ++number_finite_bounds;
    }
    if (model_.hasUb(i)) {
      mu_ += xu_[i] * zu_[i];
      ++number_finite_bounds;
    }
  }
  mu_ /= number_finite_bounds;
}
void IpmIterate::scaling() {
  scaling_.assign(model_.n_, 0.0);

  for (int i = 0; i < model_.n_; ++i) {
    if (model_.hasLb(i)) scaling_[i] += zl_[i] / xl_[i];
    if (model_.hasUb(i)) scaling_[i] += zu_[i] / xu_[i];

    // slow down the growth of theta
    if (scaling_[i] < 1e-12) scaling_[i] = sqrt(1e-12 * scaling_[i]);
  }

  // compute min and max entry in Theta
  double& min_theta = DataCollector::get()->back().min_theta;
  double& max_theta = DataCollector::get()->back().max_theta;
  min_theta = kInf;
  max_theta = 0.0;
  for (int i = 0; i < model_.n_; ++i) {
    if (scaling_[i] != 0.0) {
      min_theta = std::min(min_theta, 1.0 / scaling_[i]);
      max_theta = std::max(max_theta, 1.0 / scaling_[i]);
    }
  }
}
void IpmIterate::products() {
  double min_prod = std::numeric_limits<double>::max();
  double max_prod = 0.0;
  int num_small = 0;
  int num_large = 0;

  for (int i = 0; i < model_.n_; ++i) {
    if (model_.hasLb(i)) {
      double prod = xl_[i] * zl_[i] / mu_;
      min_prod = std::min(min_prod, prod);
      max_prod = std::max(max_prod, prod);
      if (prod < kSmallProduct) ++num_small;
      if (prod > kLargeProduct) ++num_large;
    }
    if (model_.hasUb(i)) {
      double prod = xu_[i] * zu_[i] / mu_;
      min_prod = std::min(min_prod, prod);
      max_prod = std::max(max_prod, prod);
      if (prod < kSmallProduct) ++num_small;
      if (prod > kLargeProduct) ++num_large;
    }
  }

  DataCollector::get()->back().min_prod = min_prod;
  DataCollector::get()->back().max_prod = max_prod;
  DataCollector::get()->back().num_small_prod = num_small;
  DataCollector::get()->back().num_large_prod = num_large;
}

void IpmIterate::indicators() {
  primalObj();
  dualObj();
  primalInfeas();
  dualInfeas();
  pdGap();
  products();
}

void IpmIterate::primalObj() {
  // compute unscaled primal objective
  pobj_ = dotProd(x_, model_.c_);
  if (model_.colexp_.size() > 0)
    pobj_ = std::ldexp(pobj_, -model_.cexp_ - model_.bexp_);
}
void IpmIterate::dualObj() {
  // compute unscaled dual objective
  dobj_ = dotProd(y_, model_.b_);
  for (int i = 0; i < model_.n_; ++i) {
    if (model_.hasLb(i)) dobj_ += model_.lower_[i] * zl_[i];
    if (model_.hasUb(i)) dobj_ -= model_.upper_[i] * zu_[i];
  }
  if (model_.colexp_.size() > 0)
    dobj_ = std::ldexp(dobj_, -model_.cexp_ - model_.bexp_);
}
void IpmIterate::pdGap() {
  // relative primal-dual gap
  pdgap_ = std::abs(pobj_ - dobj_) / (1 + 0.5 * std::abs(pobj_ + dobj_));
}

void IpmIterate::primalInfeas() {
  // relative infinity norm of scaled primal residuals
  pinf_ = infNorm(res1_);
  pinf_ = std::max(pinf_, infNorm(res2_));
  pinf_ = std::max(pinf_, infNorm(res3_));
  pinf_ /= (1 + model_.normScaledRhs());
}
void IpmIterate::dualInfeas() {
  // relative infinity norm of scaled dual residual
  dinf_ = infNorm(res4_) / (1 + model_.normScaledObj());
}
void IpmIterate::primalInfeasUnscaled() {
  // relative infinity norm of unscaled primal residuals
  pinf_ = 0.0;
  for (int i = 0; i < model_.m_; ++i) {
    double val = std::abs(res1_[i]);
    if (!model_.rowexp_.empty())
      val = std::ldexp(val, -model_.bexp_ - model_.rowexp_[i]);
    pinf_ = std::max(pinf_, val);
  }
  for (int i = 0; i < model_.n_; ++i) {
    double val = std::abs(res2_[i]);
    if (!model_.colexp_.empty())
      val = std::ldexp(val, -model_.bexp_ + model_.colexp_[i]);
    pinf_ = std::max(pinf_, val);

    val = std::abs(res3_[i]);
    if (!model_.colexp_.empty())
      val = std::ldexp(val, -model_.bexp_ + model_.colexp_[i]);
    pinf_ = std::max(pinf_, val);
  }
  pinf_ /= (1.0 + model_.normUnscaledRhs());
}
void IpmIterate::dualInfeasUnscaled() {
  // relative infinity norm of unscaled dual residual
  dinf_ = 0.0;
  for (int i = 0; i < model_.n_; ++i) {
    double val = std::abs(res4_[i]);
    if (model_.colexp_.size() > 0)
      val = std::ldexp(val, -model_.cexp_ - model_.colexp_[i]);
    dinf_ = std::max(dinf_, val);
  }
  dinf_ /= (1.0 + model_.normUnscaledObj());
}

void IpmIterate::residual1234() {
  // res1
  res1_ = model_.b_;
  model_.A_.alphaProductPlusY(-1.0, x_, res1_);

  // res2
  for (int i = 0; i < model_.n_; ++i) {
    if (model_.hasLb(i))
      res2_[i] = model_.lower_[i] - x_[i] + xl_[i];
    else
      res2_[i] = 0.0;
  }

  // res3
  for (int i = 0; i < model_.n_; ++i) {
    if (model_.hasUb(i))
      res3_[i] = model_.upper_[i] - x_[i] - xu_[i];
    else
      res3_[i] = 0.0;
  }

  // res4
  res4_ = model_.c_;
  model_.A_.alphaProductPlusY(-1.0, y_, res4_, true);
  for (int i = 0; i < model_.n_; ++i) {
    if (model_.hasLb(i)) res4_[i] -= zl_[i];
    if (model_.hasUb(i)) res4_[i] += zu_[i];
  }
}
void IpmIterate::residual56(double sigma) {
  for (int i = 0; i < model_.n_; ++i) {
    // res5
    if (model_.hasLb(i))
      res5_[i] = sigma * mu_ - xl_[i] * zl_[i];
    else
      res5_[i] = 0.0;

    // res6
    if (model_.hasUb(i))
      res6_[i] = sigma * mu_ - xu_[i] * zu_[i];
    else
      res6_[i] = 0.0;
  }
}

std::vector<double> IpmIterate::residual7() const {
  std::vector<double> res7(res4_);
  for (int i = 0; i < model_.n_; ++i) {
    if (model_.hasLb(i)) res7[i] -= ((res5_[i] + zl_[i] * res2_[i]) / xl_[i]);
    if (model_.hasUb(i)) res7[i] += ((res6_[i] - zu_[i] * res3_[i]) / xu_[i]);
  }
  return res7;
}
std::vector<double> IpmIterate::residual8(
    const std::vector<double>& res7) const {
  std::vector<double> res8(res1_);
  std::vector<double> temp(res7);

  // temp = (Theta^-1+Rp)^-1 * res7
  for (int i = 0; i < model_.n_; ++i)
    temp[i] /= scaling_[i] + kPrimalStaticRegularization;

  // res8 += A * temp
  model_.A_.alphaProductPlusY(1.0, temp, res8);

  return res8;
}
void IpmIterate::clearRes() {
  res1_.assign(model_.m_, 0.0);
  res2_.assign(model_.n_, 0.0);
  res3_.assign(model_.n_, 0.0);
  res4_.assign(model_.n_, 0.0);
  res5_.assign(model_.n_, 0.0);
  res6_.assign(model_.n_, 0.0);
}

void IpmIterate::prepareForUser(std::vector<double>& x, std::vector<double>& xl,
                                std::vector<double>& xu,
                                std::vector<double>& slack,
                                std::vector<double>& y, std::vector<double>& zl,
                                std::vector<double>& zu) const {
  // copy x, xl, xu, zl, zu without slacks
  x = std::vector<double>(x_.begin(), x_.begin() + model_.num_var_);
  xl = std::vector<double>(xl_.begin(), xl_.begin() + model_.num_var_);
  xu = std::vector<double>(xu_.begin(), xu_.begin() + model_.num_var_);
  zl = std::vector<double>(zl_.begin(), zl_.begin() + model_.num_var_);
  zu = std::vector<double>(zu_.begin(), zu_.begin() + model_.num_var_);

  // for the Lagrange multipliers, use slacks from zl and zu, to get correct
  // sign
  y.resize(model_.m_);
  int slack_pos = 0;
  for (int i = 0; i < model_.m_; ++i) {
    switch (model_.constraints_[i]) {
      case kConstraintTypeEqual:
        y[i] = y_[i];
        break;
      case kConstraintTypeLower:
        y[i] = zu_[model_.num_var_ + slack_pos];
        ++slack_pos;
        break;
      case kConstraintTypeUpper:
        y[i] = -zl_[model_.num_var_ + slack_pos];
        ++slack_pos;
        break;
    }
  }

  // for x-slacks, use slacks from xl and xu, to get correct sign
  slack.resize(model_.m_);
  slack_pos = 0;
  for (int i = 0; i < model_.m_; ++i) {
    switch (model_.constraints_[i]) {
      case kConstraintTypeEqual:
        slack[i] = 0.0;
        break;
      case kConstraintTypeLower:
        slack[i] = -xu_[model_.num_var_ + slack_pos];
        ++slack_pos;
        break;
      case kConstraintTypeUpper:
        slack[i] = xl_[model_.num_var_ + slack_pos];
        ++slack_pos;
        break;
    }
  }
}