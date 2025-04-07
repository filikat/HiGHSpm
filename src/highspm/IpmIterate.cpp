#include "IpmIterate.h"

#include "factorhighs/DataCollector.h"
#include "Ipm_const.h"

NewtonDir::NewtonDir(int m, int n)
    : x(n, 0.0), y(m, 0.0), xl(n, 0.0), xu(n, 0.0), zl(n, 0.0), zu(n, 0.0) {}

IpmIterate::IpmIterate(const IpmModel& model_input)
    : model{&model_input}, delta(model->m(), model->n()) {
  clearIter();
  clearRes();
}

bool IpmIterate::isNan() const {
  if (isNanVector(x) || isNanVector(xl) || isNanVector(xu) || isNanVector(y) ||
      isNanVector(zl) || isNanVector(zu))
    return true;
  return false;
}
bool IpmIterate::isInf() const {
  if (isInfVector(x) || isInfVector(xl) || isInfVector(xu) || isInfVector(y) ||
      isInfVector(zl) || isInfVector(zu))
    return true;
  return false;
}
bool IpmIterate::isResNan() const {
  if (isNanVector(res1) || isNanVector(res2) || isNanVector(res3) ||
      isNanVector(res4) || isNanVector(res5) || isNanVector(res6))
    return true;
  return false;
}
bool IpmIterate::isResInf() const {
  if (isInfVector(res1) || isInfVector(res2) || isInfVector(res3) ||
      isInfVector(res4) || isInfVector(res5) || isInfVector(res6))
    return true;
  return false;
}
bool IpmIterate::isDirNan() const {
  if (isNanVector(delta.x) || isNanVector(delta.xl) || isNanVector(delta.xu) ||
      isNanVector(delta.y) || isNanVector(delta.zl) || isNanVector(delta.zu))
    return true;
  return false;
}
bool IpmIterate::isDirInf() const {
  if (isInfVector(delta.x) || isInfVector(delta.xl) || isInfVector(delta.xu) ||
      isInfVector(delta.y) || isInfVector(delta.zl) || isInfVector(delta.zu))
    return true;
  return false;
}

void IpmIterate::computeMu() {
  mu = 0.0;
  int number_finite_bounds{};
  for (int i = 0; i < model->n(); ++i) {
    if (model->hasLb(i)) {
      mu += xl[i] * zl[i];
      ++number_finite_bounds;
    }
    if (model->hasUb(i)) {
      mu += xu[i] * zu[i];
      ++number_finite_bounds;
    }
  }
  mu /= number_finite_bounds;
}
void IpmIterate::computeScaling() {
  scaling.assign(model->n(), 0.0);

  for (int i = 0; i < model->n(); ++i) {
    if (model->hasLb(i)) scaling[i] += zl[i] / xl[i];
    if (model->hasUb(i)) scaling[i] += zu[i] / xu[i];

    // slow down the growth of theta
    if (scaling[i] < 1e-12) scaling[i] = sqrt(1e-12 * scaling[i]);
  }

  // compute min and max entry in Theta
  double& min_theta = DataCollector::get()->back().min_theta;
  double& max_theta = DataCollector::get()->back().max_theta;
  min_theta = kHighsInf;
  max_theta = 0.0;
  for (int i = 0; i < model->n(); ++i) {
    if (scaling[i] != 0.0) {
      min_theta = std::min(min_theta, 1.0 / scaling[i]);
      max_theta = std::max(max_theta, 1.0 / scaling[i]);
    }
  }
}
void IpmIterate::products() {
  double min_prod = std::numeric_limits<double>::max();
  double max_prod = 0.0;
  int num_small = 0;
  int num_large = 0;

  for (int i = 0; i < model->n(); ++i) {
    if (model->hasLb(i)) {
      double prod = xl[i] * zl[i] / mu;
      min_prod = std::min(min_prod, prod);
      max_prod = std::max(max_prod, prod);
      if (prod < kSmallProduct) ++num_small;
      if (prod > kLargeProduct) ++num_large;
    }
    if (model->hasUb(i)) {
      double prod = xu[i] * zu[i] / mu;
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
  pobj = model->offset() + dotProd(x, model->c());
}
void IpmIterate::dualObj() {
  dobj = model->offset() + dotProd(y, model->b());
  for (int i = 0; i < model->n(); ++i) {
    if (model->hasLb(i)) dobj += model->lb(i) * zl[i];
    if (model->hasUb(i)) dobj -= model->ub(i) * zu[i];
  }
}
void IpmIterate::pdGap() {
  // relative primal-dual gap
  pdgap = std::abs(pobj - dobj) / (1 + 0.5 * std::abs(pobj + dobj));
}

void IpmIterate::primalInfeas() {
  // relative infinity norm of scaled primal residuals
  pinf = infNorm(res1);
  pinf = std::max(pinf, infNorm(res2));
  pinf = std::max(pinf, infNorm(res3));
  pinf /= (1 + model->normScaledRhs());
}
void IpmIterate::dualInfeas() {
  // relative infinity norm of scaled dual residual
  dinf = infNorm(res4) / (1 + model->normScaledObj());
}
void IpmIterate::primalInfeasUnscaled() {
  // relative infinity norm of unscaled primal residuals
  pinf = 0.0;
  for (int i = 0; i < model->m(); ++i) {
    double val = std::abs(res1[i]);
    if (model->scaled()) val /= model->rowScale(i);
    pinf = std::max(pinf, val);
  }
  for (int i = 0; i < model->n(); ++i) {
    double val = std::abs(res2[i]);
    if (model->scaled()) val *= model->colScale(i);
    pinf = std::max(pinf, val);

    val = std::abs(res3[i]);
    if (model->scaled()) val *= model->colScale(i);
    pinf = std::max(pinf, val);
  }
  pinf /= (1.0 + model->normUnscaledRhs());
}
void IpmIterate::dualInfeasUnscaled() {
  // relative infinity norm of unscaled dual residual
  dinf = 0.0;
  for (int i = 0; i < model->n(); ++i) {
    double val = std::abs(res4[i]);
    if (model->scaled()) val /= model->colScale(i);
    dinf = std::max(dinf, val);
  }
  dinf /= (1.0 + model->normUnscaledObj());
}

void IpmIterate::residual1234() {
  // res1
  res1 = model->b();
  model->A().alphaProductPlusY(-1.0, x, res1);

  // res2
  for (int i = 0; i < model->n(); ++i) {
    if (model->hasLb(i))
      res2[i] = model->lb(i) - x[i] + xl[i];
    else
      res2[i] = 0.0;
  }

  // res3
  for (int i = 0; i < model->n(); ++i) {
    if (model->hasUb(i))
      res3[i] = model->ub(i) - x[i] - xu[i];
    else
      res3[i] = 0.0;
  }

  // res4
  res4 = model->c();
  model->A().alphaProductPlusY(-1.0, y, res4, true);
  for (int i = 0; i < model->n(); ++i) {
    if (model->hasLb(i)) res4[i] -= zl[i];
    if (model->hasUb(i)) res4[i] += zu[i];
  }
}
void IpmIterate::residual56(double sigma) {
  for (int i = 0; i < model->n(); ++i) {
    // res5
    if (model->hasLb(i))
      res5[i] = sigma * mu - xl[i] * zl[i];
    else
      res5[i] = 0.0;

    // res6
    if (model->hasUb(i))
      res6[i] = sigma * mu - xu[i] * zu[i];
    else
      res6[i] = 0.0;
  }
}

std::vector<double> IpmIterate::residual7() const {
  std::vector<double> res7(res4);
  for (int i = 0; i < model->n(); ++i) {
    if (model->hasLb(i)) res7[i] -= ((res5[i] + zl[i] * res2[i]) / xl[i]);
    if (model->hasUb(i)) res7[i] += ((res6[i] - zu[i] * res3[i]) / xu[i]);
  }
  return res7;
}
std::vector<double> IpmIterate::residual8(
    const std::vector<double>& res7) const {
  std::vector<double> res8(res1);
  std::vector<double> temp(res7);

  // temp = (Theta^-1+Rp)^-1 * res7
  for (int i = 0; i < model->n(); ++i)
    temp[i] /= scaling[i] + kPrimalStaticRegularization;

  // res8 += A * temp
  model->A().alphaProductPlusY(1.0, temp, res8);

  return res8;
}

void IpmIterate::clearIter() {
  x.assign(model->n(), 0.0);
  xl.assign(model->n(), 0.0);
  xu.assign(model->n(), 0.0);
  y.assign(model->m(), 0.0);
  zl.assign(model->n(), 0.0);
  zu.assign(model->n(), 0.0);
}
void IpmIterate::clearRes() {
  res1.assign(model->m(), 0.0);
  res2.assign(model->n(), 0.0);
  res3.assign(model->n(), 0.0);
  res4.assign(model->n(), 0.0);
  res5.assign(model->n(), 0.0);
  res6.assign(model->n(), 0.0);
}
void IpmIterate::clearDir() {
  delta.x.assign(model->n(), 0.0);
  delta.xl.assign(model->n(), 0.0);
  delta.xu.assign(model->n(), 0.0);
  delta.y.assign(model->m(), 0.0);
  delta.zl.assign(model->n(), 0.0);
  delta.zu.assign(model->n(), 0.0);
}

void IpmIterate::extract(std::vector<double>& x_user,
                         std::vector<double>& xl_user,
                         std::vector<double>& xu_user,
                         std::vector<double>& slack_user,
                         std::vector<double>& y_user,
                         std::vector<double>& zl_user,
                         std::vector<double>& zu_user) const {
  // Extract solution with internal format

  // Copy x, xl, xu, zl, zu without slacks
  x_user = std::vector<double>(x.begin(), x.begin() + model->n_orig());
  xl_user = std::vector<double>(xl.begin(), xl.begin() + model->n_orig());
  xu_user = std::vector<double>(xu.begin(), xu.begin() + model->n_orig());
  zl_user = std::vector<double>(zl.begin(), zl.begin() + model->n_orig());
  zu_user = std::vector<double>(zu.begin(), zu.begin() + model->n_orig());

  // For the Lagrange multipliers, use slacks from zl and zu, to get correct
  // sign. NB: there is no explicit slack stored for equality constraints.
  y_user.resize(model->m());
  int slack_pos = 0;
  for (int i = 0; i < model->m(); ++i) {
    switch (model->constraint(i)) {
      case '=':
        y_user[i] = y[i];
        break;
      case '>':
        y_user[i] = zu[model->n_orig() + slack_pos];
        ++slack_pos;
        break;
      case '<':
        y_user[i] = -zl[model->n_orig() + slack_pos];
        ++slack_pos;
        break;
    }
  }

  // For x-slacks, use slacks from xl and xu, to get correct sign.
  // NB: there is no explicit slack stored for equality constraints.
  slack_user.resize(model->m());
  slack_pos = 0;
  for (int i = 0; i < model->m(); ++i) {
    switch (model->constraint(i)) {
      case '=':
        slack_user[i] = 0.0;
        break;
      case '>':
        slack_user[i] = -xu[model->n_orig() + slack_pos];
        ++slack_pos;
        break;
      case '<':
        slack_user[i] = xl[model->n_orig() + slack_pos];
        ++slack_pos;
        break;
    }
  }
}

void IpmIterate::extract(std::vector<double>& x_user,
                         std::vector<double>& slack_user,
                         std::vector<double>& y_user,
                         std::vector<double>& z_user) const {
  // Extract solution with format for crossover

  // Construct complementary point (x_temp, y_temp, z_temp)
  std::vector<double> x_temp, y_temp, z_temp;
  dropToComplementarity(x_temp, y_temp, z_temp);

  // Both x_temp and z_temp include slacks.
  // They are removed from x and z, but they are used to compute slack and y.

  // Remove slacks from x and z
  x_user =
      std::vector<double>(x_temp.begin(), x_temp.begin() + model->n_orig());
  z_user =
      std::vector<double>(z_temp.begin(), z_temp.begin() + model->n_orig());

  // For inequality constraints, the corresponding z-slack may have been dropped
  // to zero, so build y from z-slacks.
  // NB: there is no explicit slack stored for equality constraints.
  y_user.resize(model->m());
  int slack_pos = 0;
  for (int i = 0; i < model->m(); ++i) {
    switch (model->constraint(i)) {
      case '=':
        y_user[i] = y_temp[i];
        break;
      case '>':
      case '<':
        y_user[i] = -z_temp[model->n_orig() + slack_pos];
        ++slack_pos;
        break;
    }
  }

  // Use slacks from x_temp and add slack for equality constraints.
  // NB: there is no explicit slack stored for equality constraints.
  slack_user.resize(model->m());
  slack_pos = 0;
  for (int i = 0; i < model->m(); ++i) {
    switch (model->constraint(i)) {
      case '=':
        slack_user[i] = 0.0;
        break;
      case '>':
      case '<':
        slack_user[i] = x_temp[model->n_orig() + slack_pos];
        ++slack_pos;
        break;
    }
  }
}

void IpmIterate::dropToComplementarity(std::vector<double>& x_cmp,
                                       std::vector<double>& y_cmp,
                                       std::vector<double>& z_cmp) const {
  x_cmp.assign(model->n(), 0.0);
  z_cmp.assign(model->n(), 0.0);
  y_cmp = y;

  for (int j = 0; j < model->n(); ++j) {
    // value of x_[j] within bounds
    double xj = std::max(x[j], model->lb(j));
    xj = std::min(xj, model->ub(j));

    // FIXED VARIABLE
    if (model->lb(j) == model->ub(j)) {
      x_cmp[j] = model->lb(j);
      z_cmp[j] = zl[j] - zu[j];
    }

    // BOTH BOUNDS FINITE
    else if (model->hasLb(j) && model->hasUb(j)) {
      if (zl[j] * xu[j] >= zu[j] * xl[j]) {
        // xlj/zlj <= xuj/zuj
        // Primal lower is smaller than primal upper, wrt respective duals
        if (zl[j] >= xl[j]) {
          // drop x to lower bound, set z positive
          x_cmp[j] = model->lb(j);
          z_cmp[j] = std::max(0.0, zl[j] - zu[j]);
        } else {
          // drop z to zero, set x within bounds
          x_cmp[j] = xj;
          z_cmp[j] = 0.0;
        }
      } else {
        // xuj/zuj < xlj/zlj
        // Primal upper is smaller than primal lower, wrt respective duals
        if (zu[j] >= xu[j]) {
          // drop x to upper bound, set z negative
          x_cmp[j] = model->ub(j);
          z_cmp[j] = std::min(0.0, zl[j] - zu[j]);
        } else {
          // drop z to zero, set x within bounds
          x_cmp[j] = xj;
          z_cmp[j] = 0.0;
        }
      }
    }

    // LOWER BOUND FINITE
    else if (model->hasLb(j)) {
      if (zl[j] >= xl[j]) {
        // drop x to lower bound, set z positive
        x_cmp[j] = model->lb(j);
        z_cmp[j] = std::max(0.0, zl[j] - zu[j]);
      } else {
        // drop z to zero, set x within bounds
        x_cmp[j] = xj;
        z_cmp[j] = 0.0;
      }
    }

    // UPPER BOUND FINITE
    else if (model->hasUb(j)) {
      if (zu[j] >= xu[j]) {
        // drop x to upper bound, set z negative
        x_cmp[j] = model->ub(j);
        z_cmp[j] = std::min(0.0, zl[j] - zu[j]);
      } else {
        // drop z to zero, set x within bounds
        x_cmp[j] = xj;
        z_cmp[j] = 0.0;
      }
    }

    // NO BOUNDS
    else {
      x_cmp[j] = xj;
      z_cmp[j] = 0.0;
    }
  }
}

double IpmIterate::infeasAfterDropping() const {
  // Compute estimate of residuals after dropping to complementarity (taken from
  // ipx).

  double pinf_max = 0.0;
  double dinf_max = 0.0;

  for (int j = 0; j < model->n(); ++j) {
    double xdrop = 0.0;
    double zdrop = 0.0;

    if (model->hasLb(j) && model->hasUb(j)) {
      // BOTH BOUNDS FINITE
      if (zl[j] * xu[j] >= zu[j] * xl[j]) {
        if (zl[j] >= xl[j])
          xdrop = x[j] - model->lb(j);
        else
          zdrop = zl[j] - zu[j];
      } else {
        if (zu[j] >= xu[j])
          xdrop = x[j] - model->ub(j);
        else
          zdrop = zl[j] - zu[j];
      }
    }

    // LOWER BOUND FINITE
    else if (model->hasLb(j)) {
      if (zl[j] >= xl[j])
        xdrop = x[j] - model->lb(j);
      else
        zdrop = zl[j] - zu[j];
    }

    // UPPER BOUND FINITE
    else if (model->hasUb(j)) {
      if (zu[j] >= xu[j])
        xdrop = x[j] - model->ub(j);
      else
        zdrop = zl[j] - zu[j];
    }

    // largest entry in column j of A
    double Amax = 0.0;
    for (int el = model->A().start_[j]; el < model->A().start_[j + 1]; ++el)
      Amax = std::max(Amax, std::abs(model->A().value_[el]));

    pinf_max = std::max(pinf_max, std::abs(xdrop) * Amax);
    dinf_max = std::max(dinf_max, std::abs(zdrop));
  }

  pinf_max /= (1.0 + model->normScaledRhs());
  dinf_max /= (1.0 + model->normScaledObj());

  return std::max(pinf_max, dinf_max);
}