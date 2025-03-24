#include "IpmIterate.h"

#include <cassert>

#include "../FactorHiGHS/DataCollector.h"
#include "Ipm_const.h"

NewtonDir::NewtonDir(int m, int n)
    : x(n, 0.0), y(m, 0.0), xl(n, 0.0), xu(n, 0.0), zl(n, 0.0), zu(n, 0.0) {}

IpmIterate::IpmIterate(const IpmModel& model)
    : model_{model},
      delta_(model_.m(), model_.n()),
      dx_{delta_.x},
      dxl_{delta_.xl},
      dxu_{delta_.xu},
      dy_{delta_.y},
      dzl_{delta_.zl},
      dzu_{delta_.zu} {
  clearIter();
  clearRes();
}

bool IpmIterate::isNan() const {
  if (isNanVector(x_) || isNanVector(xl_) || isNanVector(xu_) ||
      isNanVector(y_) || isNanVector(zl_) || isNanVector(zu_))
    return true;
  return false;
}
bool IpmIterate::isInf() const {
  if (isInfVector(x_) || isInfVector(xl_) || isInfVector(xu_) ||
      isInfVector(y_) || isInfVector(zl_) || isInfVector(zu_))
    return true;
  return false;
}
bool IpmIterate::isResNan() const {
  if (isNanVector(res1_) || isNanVector(res2_) || isNanVector(res3_) ||
      isNanVector(res4_) || isNanVector(res5_) || isNanVector(res6_))
    return true;
  return false;
}
bool IpmIterate::isResInf() const {
  if (isInfVector(res1_) || isInfVector(res2_) || isInfVector(res3_) ||
      isInfVector(res4_) || isInfVector(res5_) || isInfVector(res6_))
    return true;
  return false;
}
bool IpmIterate::isDirNan() const {
  if (isNanVector(dx_) || isNanVector(dxl_) || isNanVector(dxu_) ||
      isNanVector(dy_) || isNanVector(dzl_) || isNanVector(dzu_))
    return true;
  return false;
}
bool IpmIterate::isDirInf() const {
  if (isInfVector(dx_) || isInfVector(dxl_) || isInfVector(dxu_) ||
      isInfVector(dy_) || isInfVector(dzl_) || isInfVector(dzu_))
    return true;
  return false;
}

void IpmIterate::mu() {
  mu_ = 0.0;
  int number_finite_bounds{};
  for (int i = 0; i < model_.n(); ++i) {
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
  scaling_.assign(model_.n(), 0.0);

  for (int i = 0; i < model_.n(); ++i) {
    if (model_.hasLb(i)) scaling_[i] += zl_[i] / xl_[i];
    if (model_.hasUb(i)) scaling_[i] += zu_[i] / xu_[i];

    // slow down the growth of theta
    if (scaling_[i] < 1e-12) scaling_[i] = sqrt(1e-12 * scaling_[i]);
  }

  // compute min and max entry in Theta
  double& min_theta = DataCollector::get()->back().min_theta;
  double& max_theta = DataCollector::get()->back().max_theta;
  min_theta = kHighsInf;
  max_theta = 0.0;
  for (int i = 0; i < model_.n(); ++i) {
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

  for (int i = 0; i < model_.n(); ++i) {
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
  pobj_ = model_.offset() + dotProd(x_, model_.c());
}
void IpmIterate::dualObj() {
  dobj_ = model_.offset() + dotProd(y_, model_.b());
  for (int i = 0; i < model_.n(); ++i) {
    if (model_.hasLb(i)) dobj_ += model_.lb(i) * zl_[i];
    if (model_.hasUb(i)) dobj_ -= model_.ub(i) * zu_[i];
  }
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
  for (int i = 0; i < model_.m(); ++i) {
    double val = std::abs(res1_[i]);
    if (model_.scaled()) val /= model_.rowScale(i);
    pinf_ = std::max(pinf_, val);
  }
  for (int i = 0; i < model_.n(); ++i) {
    double val = std::abs(res2_[i]);
    if (model_.scaled()) val *= model_.colScale(i);
    pinf_ = std::max(pinf_, val);

    val = std::abs(res3_[i]);
    if (model_.scaled()) val *= model_.colScale(i);
    pinf_ = std::max(pinf_, val);
  }
  pinf_ /= (1.0 + model_.normUnscaledRhs());
}
void IpmIterate::dualInfeasUnscaled() {
  // relative infinity norm of unscaled dual residual
  dinf_ = 0.0;
  for (int i = 0; i < model_.n(); ++i) {
    double val = std::abs(res4_[i]);
    if (model_.scaled()) val /= model_.colScale(i);
    dinf_ = std::max(dinf_, val);
  }
  dinf_ /= (1.0 + model_.normUnscaledObj());
}

void IpmIterate::residual1234() {
  // res1
  res1_ = model_.b();
  model_.A().alphaProductPlusY(-1.0, x_, res1_);

  // res2
  for (int i = 0; i < model_.n(); ++i) {
    if (model_.hasLb(i))
      res2_[i] = model_.lb(i) - x_[i] + xl_[i];
    else
      res2_[i] = 0.0;
  }

  // res3
  for (int i = 0; i < model_.n(); ++i) {
    if (model_.hasUb(i))
      res3_[i] = model_.ub(i) - x_[i] - xu_[i];
    else
      res3_[i] = 0.0;
  }

  // res4
  res4_ = model_.c();
  model_.A().alphaProductPlusY(-1.0, y_, res4_, true);
  for (int i = 0; i < model_.n(); ++i) {
    if (model_.hasLb(i)) res4_[i] -= zl_[i];
    if (model_.hasUb(i)) res4_[i] += zu_[i];
  }

  assert(!isResNan() && !isResInf());
}
void IpmIterate::residual56(double sigma) {
  for (int i = 0; i < model_.n(); ++i) {
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

  assert(!isResNan() && !isResInf());
}

std::vector<double> IpmIterate::residual7() const {
  std::vector<double> res7(res4_);
  for (int i = 0; i < model_.n(); ++i) {
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
  for (int i = 0; i < model_.n(); ++i)
    temp[i] /= scaling_[i] + kPrimalStaticRegularization;

  // res8 += A * temp
  model_.A().alphaProductPlusY(1.0, temp, res8);

  return res8;
}

void IpmIterate::clearIter() {
  x_.assign(model_.n(), 0.0);
  xl_.assign(model_.n(), 0.0);
  xu_.assign(model_.n(), 0.0);
  y_.assign(model_.m(), 0.0);
  zl_.assign(model_.n(), 0.0);
  zu_.assign(model_.n(), 0.0);
}
void IpmIterate::clearRes() {
  res1_.assign(model_.m(), 0.0);
  res2_.assign(model_.n(), 0.0);
  res3_.assign(model_.n(), 0.0);
  res4_.assign(model_.n(), 0.0);
  res5_.assign(model_.n(), 0.0);
  res6_.assign(model_.n(), 0.0);
}
void IpmIterate::clearDir() {
  dx_.assign(model_.n(), 0.0);
  dxl_.assign(model_.n(), 0.0);
  dxu_.assign(model_.n(), 0.0);
  dy_.assign(model_.m(), 0.0);
  dzl_.assign(model_.n(), 0.0);
  dzu_.assign(model_.n(), 0.0);
}

void IpmIterate::extract(std::vector<double>& x, std::vector<double>& xl,
                         std::vector<double>& xu, std::vector<double>& slack,
                         std::vector<double>& y, std::vector<double>& zl,
                         std::vector<double>& zu) const {
  // Extract solution with internal format

  // Copy x, xl, xu, zl, zu without slacks
  x = std::vector<double>(x_.begin(), x_.begin() + model_.n_orig());
  xl = std::vector<double>(xl_.begin(), xl_.begin() + model_.n_orig());
  xu = std::vector<double>(xu_.begin(), xu_.begin() + model_.n_orig());
  zl = std::vector<double>(zl_.begin(), zl_.begin() + model_.n_orig());
  zu = std::vector<double>(zu_.begin(), zu_.begin() + model_.n_orig());

  // For the Lagrange multipliers, use slacks from zl and zu, to get correct
  // sign. NB: there is no explicit slack stored for equality constraints.
  y.resize(model_.m());
  int slack_pos = 0;
  for (int i = 0; i < model_.m(); ++i) {
    switch (model_.constraint(i)) {
      case '=':
        y[i] = y_[i];
        break;
      case '>':
        y[i] = zu_[model_.n_orig() + slack_pos];
        ++slack_pos;
        break;
      case '<':
        y[i] = -zl_[model_.n_orig() + slack_pos];
        ++slack_pos;
        break;
    }
  }

  // For x-slacks, use slacks from xl and xu, to get correct sign.
  // NB: there is no explicit slack stored for equality constraints.
  slack.resize(model_.m());
  slack_pos = 0;
  for (int i = 0; i < model_.m(); ++i) {
    switch (model_.constraint(i)) {
      case '=':
        slack[i] = 0.0;
        break;
      case '>':
        slack[i] = -xu_[model_.n_orig() + slack_pos];
        ++slack_pos;
        break;
      case '<':
        slack[i] = xl_[model_.n_orig() + slack_pos];
        ++slack_pos;
        break;
    }
  }
}

void IpmIterate::extract(std::vector<double>& x, std::vector<double>& slack,
                         std::vector<double>& y, std::vector<double>& z) const {
  // Extract solution with format for crossover

  // Construct complementary point (x_temp, y_temp, z_temp)
  std::vector<double> x_temp, y_temp, z_temp;
  dropToComplementarity(x_temp, y_temp, z_temp);

  // Both x_temp and z_temp include slacks.
  // They are removed from x and z, but they are used to compute slack and y.

  // Remove slacks from x and z
  x = std::vector<double>(x_temp.begin(), x_temp.begin() + model_.n_orig());
  z = std::vector<double>(z_temp.begin(), z_temp.begin() + model_.n_orig());

  // For inequality constraints, the corresponding z-slack may have been dropped
  // to zero, so build y from z-slacks.
  // NB: there is no explicit slack stored for equality constraints.
  y.resize(model_.m());
  int slack_pos = 0;
  for (int i = 0; i < model_.m(); ++i) {
    switch (model_.constraint(i)) {
      case '=':
        y[i] = y_temp[i];
        break;
      case '>':
      case '<':
        y[i] = -z_temp[model_.n_orig() + slack_pos];
        ++slack_pos;
        break;
    }
  }

  // Use slacks from x_temp and add slack for equality constraints.
  // NB: there is no explicit slack stored for equality constraints.
  slack.resize(model_.m());
  slack_pos = 0;
  for (int i = 0; i < model_.m(); ++i) {
    switch (model_.constraint(i)) {
      case '=':
        slack[i] = 0.0;
        break;
      case '>':
      case '<':
        slack[i] = x_temp[model_.n_orig() + slack_pos];
        ++slack_pos;
        break;
    }
  }
}

void IpmIterate::dropToComplementarity(std::vector<double>& x,
                                       std::vector<double>& y,
                                       std::vector<double>& z) const {
  x.assign(model_.n(), 0.0);
  z.assign(model_.n(), 0.0);
  y = y_;

  for (int j = 0; j < model_.n(); ++j) {
    // value of x_[j] within bounds
    double xj = std::max(x_[j], model_.lb(j));
    xj = std::min(xj, model_.ub(j));

    // FIXED VARIABLE
    if (model_.lb(j) == model_.ub(j)) {
      x[j] = model_.lb(j);
      z[j] = zl_[j] - zu_[j];
    }

    // BOTH BOUNDS FINITE
    else if (model_.hasLb(j) && model_.hasUb(j)) {
      if (zl_[j] * xu_[j] >= zu_[j] * xl_[j]) {
        // xlj/zlj <= xuj/zuj
        // Primal lower is smaller than primal upper, wrt respective duals
        if (zl_[j] >= xl_[j]) {
          // drop x to lower bound, set z positive
          x[j] = model_.lb(j);
          z[j] = std::max(0.0, zl_[j] - zu_[j]);
        } else {
          // drop z to zero, set x within bounds
          x[j] = xj;
          z[j] = 0.0;
        }
      } else {
        // xuj/zuj < xlj/zlj
        // Primal upper is smaller than primal lower, wrt respective duals
        if (zu_[j] >= xu_[j]) {
          // drop x to upper bound, set z negative
          x[j] = model_.ub(j);
          z[j] = std::min(0.0, zl_[j] - zu_[j]);
        } else {
          // drop z to zero, set x within bounds
          x[j] = xj;
          z[j] = 0.0;
        }
      }
    }

    // LOWER BOUND FINITE
    else if (model_.hasLb(j)) {
      if (zl_[j] >= xl_[j]) {
        // drop x to lower bound, set z positive
        x[j] = model_.lb(j);
        z[j] = std::max(0.0, zl_[j] - zu_[j]);
      } else {
        // drop z to zero, set x within bounds
        x[j] = xj;
        z[j] = 0.0;
      }
    }

    // UPPER BOUND FINITE
    else if (model_.hasUb(j)) {
      if (zu_[j] >= xu_[j]) {
        // drop x to upper bound, set z negative
        x[j] = model_.ub(j);
        z[j] = std::min(0.0, zl_[j] - zu_[j]);
      } else {
        // drop z to zero, set x within bounds
        x[j] = xj;
        z[j] = 0.0;
      }
    }

    // NO BOUNDS
    else {
      x[j] = xj;
      z[j] = 0.0;
    }
  }
}

double IpmIterate::infeasAfterDropping() const {
  // Compute estimate of residuals after dropping to complementarity (taken from
  // ipx).

  double pinf_max = 0.0;
  double dinf_max = 0.0;

  for (int j = 0; j < model_.n(); ++j) {
    double xdrop = 0.0;
    double zdrop = 0.0;

    if (model_.hasLb(j) && model_.hasUb(j)) {
      // BOTH BOUNDS FINITE
      if (zl_[j] * xu_[j] >= zu_[j] * xl_[j]) {
        if (zl_[j] >= xl_[j])
          xdrop = x_[j] - model_.lb(j);
        else
          zdrop = zl_[j] - zu_[j];
      } else {
        if (zu_[j] >= xu_[j])
          xdrop = x_[j] - model_.ub(j);
        else
          zdrop = zl_[j] - zu_[j];
      }
    }

    // LOWER BOUND FINITE
    else if (model_.hasLb(j)) {
      if (zl_[j] >= xl_[j])
        xdrop = x_[j] - model_.lb(j);
      else
        zdrop = zl_[j] - zu_[j];
    }

    // UPPER BOUND FINITE
    else if (model_.hasUb(j)) {
      if (zu_[j] >= xu_[j])
        xdrop = x_[j] - model_.ub(j);
      else
        zdrop = zl_[j] - zu_[j];
    }

    // largest entry in column j of A
    double Amax = 0.0;
    for (int el = model_.A().start_[j]; el < model_.A().start_[j + 1]; ++el)
      Amax = std::max(Amax, std::abs(model_.A().value_[el]));

    pinf_max = std::max(pinf_max, std::abs(xdrop) * Amax);
    dinf_max = std::max(dinf_max, std::abs(zdrop));
  }

  pinf_max /= (1.0 + model_.normScaledRhs());
  dinf_max /= (1.0 + model_.normScaledObj());

  return std::max(pinf_max, dinf_max);
}