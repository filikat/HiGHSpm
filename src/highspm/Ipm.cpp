#include "Ipm.h"

#include <cassert>
#include <cmath>
#include <iostream>

#include "parallel/HighsParallel.h"

namespace highspm {

void Ipm::load(const Int num_var, const Int num_con, const double* obj,
               const double* rhs, const double* lower, const double* upper,
               const Int* A_ptr, const Int* A_rows, const double* A_vals,
               const char* constraints, double offset,
               const std::string& pb_name) {
  if (!obj || !rhs || !lower || !upper || !A_ptr || !A_rows || !A_vals ||
      !constraints)
    return;

  model_.init(num_var, num_con, obj, rhs, lower, upper, A_ptr, A_rows, A_vals,
              constraints, offset, pb_name);

  m_ = model_.m();
  n_ = model_.n();

  info_.ipx_used = false;
  info_.m_solver = m_;
  info_.n_solver = n_;
  info_.m_original = num_con;
  info_.n_original = num_var;
}

void Ipm::setOptions(const Options& options) { options_ = options; }

void Ipm::solve() {
  if (!model_.ready()) {
    info_.ipm_status = kIpmStatusError;
    return;
  }

  DataCollector::initialise();
  printInfo();

  runIpm();
  refineWithIpx();

  DataCollector::get()->printTimes();
  DataCollector::get()->printIter();
  DataCollector::terminate();
}

void Ipm::runIpm() {
  if (initialise()) return;

  while (iter_ < options_.max_iter) {
    if (prepareIter()) break;
    if (predictor()) break;
    if (correctors()) break;
    makeStep();
  }

  terminate();
}

bool Ipm::initialise() {
  // Prepare ipm for execution.
  // Return true if an error occurred.

  // start timer
  clock_.start();

  // initialise iterate object
  it_.reset(new IpmIterate(model_));

  // initialise linear solver
  LS_.reset(new FactorHiGHSSolver(options_, &info_));
  if (LS_->setup(model_.A(), options_)) {
    info_.ipm_status = kIpmStatusError;
    return true;
  }
  LS_->clear();

  // decide number of correctors to use
  maxCorrectors();

  startingPoint();

  it_->residual1234();
  it_->computeMu();
  it_->indicators();

  printOutput();

  return false;
}

void Ipm::terminate() {
  info_.ipm_iter = iter_;
  if (info_.ipm_status == kIpmStatusNotRun)
    info_.ipm_status = kIpmStatusMaxIter;

  info_.option_nla = options_.nla;
  info_.option_par = options_.parallel;
}

bool Ipm::prepareIter() {
  // Prepare next iteration.
  // Return true if Ipm main loop should be stopped

  if (checkIterate()) return true;
  if (checkBadIter()) return true;
  if (checkTermination()) return true;

  ++iter_;

  // Clear Newton direction
  it_->clearDir();

  // Clear any existing data in the linear solver
  LS_->clear();

  // compute theta inverse
  it_->computeScaling();

  return false;
}

bool Ipm::predictor() {
  // Compute affine scaling direction.
  // Return true if an error occurred.

  // compute sigma and residuals for affine scaling direction
  sigmaAffine();
  it_->residual56(sigma_);

  if (solveNewtonSystem(it_->delta)) return true;
  if (recoverDirection(it_->delta)) return true;

  return false;
}

bool Ipm::correctors() {
  // Compute multiple centrality correctors.
  // Return true if an error occurred.

  sigmaCorrectors();
  if (centralityCorrectors()) return true;

  return false;
}

bool Ipm::prepareIpx() {
  // Return true if an error occurred;

  ipx::Parameters ipx_param;
  ipx_param.display = 1;
  ipx_param.dualize = 0;
  ipx_param.run_crossover = options_.crossover;
  ipx_param.ipm_feasibility_tol = options_.feasibility_tol;
  ipx_param.ipm_optimality_tol = options_.optimality_tol;
  ipx_lps_.SetParameters(ipx_param);

  Int load_status = model_.loadIntoIpx(ipx_lps_);

  if (load_status) {
    printf("Error loading model into IPX\n");
    return true;
  }

  std::vector<double> x, xl, xu, slack, y, zl, zu;
  getSolution(x, xl, xu, slack, y, zl, zu);

  Int start_point_status = ipx_lps_.LoadIPMStartingPoint(
      x.data(), xl.data(), xu.data(), slack.data(), y.data(), zl.data(),
      zu.data());

  if (start_point_status) {
    printf("Error loading starting point into IPX\n");
    return true;
  }

  return false;
}

void Ipm::refineWithIpx() {
  if (info_.ipm_status == kIpmStatusError) return;

  if (info_.ipm_status < kIpmStatusOptimal && options_.refine_with_ipx) {
    printf("\nIpm did not converge, restarting with IPX\n\n");
  } else if (options_.crossover == kOptionCrossoverOn) {
    printf("\nIpm converged, running crossover with IPX\n\n");
  } else {
    return;
  }

  if (prepareIpx()) return;

  ipx_lps_.Solve();

  info_.ipx_info = ipx_lps_.GetInfo();

  if (info_.ipx_info.status_ipm == IPX_STATUS_optimal)
    info_.ipm_status = kIpmStatusPDFeas;

  if (info_.ipx_info.status_crossover == IPX_STATUS_optimal)
    info_.ipm_status = kIpmStatusBasic;

  info_.ipx_used = true;
}

void Ipm::runCrossover() {
  prepareIpx();

  std::vector<double> x, slack, y, z;
  getSolution(x, slack, y, z);

  ipx_lps_.CrossoverFromStartingPoint(x.data(), slack.data(), y.data(),
                                      z.data());

  info_.ipx_info = ipx_lps_.GetInfo();
  info_.ipx_used = true;
}

bool Ipm::solveNewtonSystem(NewtonDir& delta) {
  std::vector<double>& theta_inv = it_->scaling;

  std::vector<double> res7 = it_->residual7();

  // NORMAL EQUATIONS
  if (options_.nla == kOptionNlaNormEq) {
    std::vector<double> res8 = it_->residual8(res7);

    // factorise normal equations, if not yet done
    if (!LS_->valid_ && LS_->factorNE(model_.A(), theta_inv)) goto failure;

    // solve with normal equations
    if (LS_->solveNE(res8, delta.y)) goto failure;

    // Compute delta.x
    // Deltax = A^T * Deltay - res7;
    delta.x = res7;
    model_.A().alphaProductPlusY(-1.0, delta.y, delta.x, true);
    vectorScale(delta.x, -1.0);

    // Deltax = (Theta^-1+Rp)^-1 * Deltax
    for (Int i = 0; i < n_; ++i)
      delta.x[i] /= theta_inv[i] + kPrimalStaticRegularisation;

  }

  // AUGMENTED SYSTEM
  else {
    // factorise augmented system, if not yet done
    if (!LS_->valid_ && LS_->factorAS(model_.A(), theta_inv)) goto failure;

    // solve with augmented system
    if (LS_->solveAS(res7, it_->res1, delta.x, delta.y)) goto failure;
  }

  return false;

// Failure occured in factorisation or solve
failure:
  std::cerr << "Error while solving Newton system\n";
  info_.ipm_status = kIpmStatusError;
  return true;
}

bool Ipm::recoverDirection(NewtonDir& delta) {
  // Recover components xl, xu, zl, zu of partial direction delta.
  std::vector<double>& xl = it_->xl;
  std::vector<double>& xu = it_->xu;
  std::vector<double>& zl = it_->zl;
  std::vector<double>& zu = it_->zu;
  std::vector<double>& res2 = it_->res2;
  std::vector<double>& res3 = it_->res3;
  std::vector<double>& res4 = it_->res4;
  std::vector<double>& res5 = it_->res5;
  std::vector<double>& res6 = it_->res6;

  for (Int i = 0; i < n_; ++i) {
    if (model_.hasLb(i) || model_.hasUb(i)) {
      delta.xl[i] = delta.x[i] - res2[i];
      delta.zl[i] = (res5[i] - zl[i] * delta.xl[i]) / xl[i];
    } else {
      delta.xl[i] = 0.0;
      delta.zl[i] = 0.0;
    }
  }
  for (Int i = 0; i < n_; ++i) {
    if (model_.hasLb(i) || model_.hasUb(i)) {
      delta.xu[i] = res3[i] - delta.x[i];
      delta.zu[i] = (res6[i] - zu[i] * delta.xu[i]) / xu[i];
    } else {
      delta.xu[i] = 0.0;
      delta.zu[i] = 0.0;
    }
  }

  // not sure if this has any effect, but IPX uses it
  std::vector<double> Atdy(n_);
  model_.A().alphaProductPlusY(1.0, delta.y, Atdy, true);
  for (Int i = 0; i < n_; ++i) {
    if (model_.hasLb(i) || model_.hasUb(i)) {
      if (std::isfinite(xl[i]) && std::isfinite(xu[i])) {
        if (zl[i] * xu[i] >= zu[i] * xl[i])
          delta.zl[i] = res4[i] + delta.zu[i] - Atdy[i];
        else
          delta.zu[i] = -res4[i] + delta.zl[i] + Atdy[i];
      } else if (std::isfinite(xl[i])) {
        delta.zl[i] = res4[i] + delta.zu[i] - Atdy[i];
      } else {
        delta.zu[i] = -res4[i] + delta.zl[i] + Atdy[i];
      }
    }
  }

  backwardError(delta);

  // Check for NaN of Inf
  if (it_->isDirNan()) {
    std::cerr << "Direction is nan\n";
    info_.ipm_status = kIpmStatusError;
    return true;
  } else if (it_->isDirInf()) {
    std::cerr << "Direction is inf\n";
    info_.ipm_status = kIpmStatusError;
    return true;
  }
  return false;
}

double Ipm::stepToBoundary(const std::vector<double>& x,
                           const std::vector<double>& dx,
                           const std::vector<double>* cor, double weight,
                           bool lo, Int* block) const {
  // Compute the largest alpha s.t. x + alpha * dx >= 0.
  // If cor is valid, consider x + alpha * (dx + w * cor) instead.
  // Use lo=1 for xl and zl, lo=0 for xu and zu.
  // Return the blocking index in block.

  const double damp = 1.0 - std::numeric_limits<double>::epsilon();

  double alpha = 1.0;
  Int bl = -1;

  for (Int i = 0; i < x.size(); ++i) {
    if ((lo && model_.hasLb(i)) || (!lo && model_.hasUb(i))) {
      double c = (cor ? (*cor)[i] * weight : 0.0);
      if (x[i] + alpha * (dx[i] + c) < 0.0) {
        alpha = -(x[i] * damp) / (dx[i] + c);
        bl = i;
      }
    }
  }
  if (block) *block = bl;
  return alpha;
}

void Ipm::stepsToBoundary(double& alpha_primal, double& alpha_dual,
                          const NewtonDir& delta, const NewtonDir* cor,
                          double weight) const {
  // compute primal and dual steps to boundary, given direction, corrector and
  // weight.

  Int block;
  double axl = stepToBoundary(it_->xl, delta.xl, cor ? &(cor->xl) : nullptr,
                              weight, true);
  double axu = stepToBoundary(it_->xu, delta.xu, cor ? &(cor->xu) : nullptr,
                              weight, false);
  double azl = stepToBoundary(it_->zl, delta.zl, cor ? &(cor->zl) : nullptr,
                              weight, true);
  double azu = stepToBoundary(it_->zu, delta.zu, cor ? &(cor->zu) : nullptr,
                              weight, false);

  alpha_primal = std::min(axl, axu);
  alpha_primal = std::min(alpha_primal, 1.0);
  alpha_dual = std::min(azl, azu);
  alpha_dual = std::min(alpha_dual, 1.0);
}

void Ipm::stepSizes() {
  // Compute primal and dual stepsizes.
  std::vector<double>& xl = it_->xl;
  std::vector<double>& xu = it_->xu;
  std::vector<double>& zl = it_->zl;
  std::vector<double>& zu = it_->zu;
  std::vector<double>& dxl = it_->delta.xl;
  std::vector<double>& dxu = it_->delta.xu;
  std::vector<double>& dzl = it_->delta.zl;
  std::vector<double>& dzu = it_->delta.zu;

  // parameters for Mehrotra heuristic
  const double gamma_f = 0.9;
  const double gamma_a = 1.0 / (1.0 - gamma_f);

  // compute stepsizes and blocking components
  Int block_xl, block_xu, block_zl, block_zu;
  double alpha_xl = stepToBoundary(xl, dxl, nullptr, 0, true, &block_xl);
  double alpha_xu = stepToBoundary(xu, dxu, nullptr, 0, false, &block_xu);
  double alpha_zl = stepToBoundary(zl, dzl, nullptr, 0, true, &block_zl);
  double alpha_zu = stepToBoundary(zu, dzu, nullptr, 0, false, &block_zu);

  double max_p = std::min(alpha_xl, alpha_xu);
  double max_d = std::min(alpha_zl, alpha_zu);

  // compute mu with current stepsizes
  double mu_full = 0.0;
  Int num_finite = 0;
  for (Int i = 0; i < n_; ++i) {
    if (model_.hasLb(i)) {
      mu_full += (xl[i] + max_p * dxl[i]) * (zl[i] + max_d * dzl[i]);
      ++num_finite;
    }
    if (model_.hasUb(i)) {
      mu_full += (xu[i] + max_p * dxu[i]) * (zu[i] + max_d * dzu[i]);
      ++num_finite;
    }
  }
  mu_full /= num_finite;
  mu_full /= gamma_a;

  // compute new stepsizes based on Mehrotra heuristic

  // primal
  double alpha_p = 1.0;
  Int block_p = -1;
  if (max_p < 1.0) {
    if (alpha_xl <= alpha_xu) {
      block_p = block_xl;
      double temp = mu_full / (zl[block_p] + max_d * dzl[block_p]);
      alpha_p = (temp - xl[block_p]) / dxl[block_p];
    } else {
      block_p = block_xu;
      double temp = mu_full / (zu[block_p] + max_d * dzu[block_p]);
      alpha_p = (temp - xu[block_p]) / dxu[block_p];
    }
    alpha_p = std::max(alpha_p, gamma_f * max_p);
    alpha_p = std::min(alpha_p, 1.0);
    assert(block_p >= 0);
  }

  // dual
  double alpha_d = 1.0;
  Int block_d = -1;
  if (max_d < 1.0) {
    if (alpha_zl <= alpha_zu) {
      block_d = block_zl;
      double temp = mu_full / (xl[block_d] + max_p * dxl[block_d]);
      alpha_d = (temp - zl[block_d]) / dzl[block_d];
    } else {
      block_d = block_zu;
      double temp = mu_full / (xu[block_d] + max_p * dxu[block_d]);
      alpha_d = (temp - zu[block_d]) / dzu[block_d];
    }
    alpha_d = std::max(alpha_d, gamma_f * max_d);
    alpha_d = std::min(alpha_d, 1.0);
    assert(block_d >= 0);
  }

  alpha_primal_ = std::min(alpha_p, 1.0 - 1e-4);
  alpha_dual_ = std::min(alpha_d, 1.0 - 1e-4);

  assert(alpha_primal_ > 0 && alpha_primal_ < 1 && alpha_dual_ > 0 &&
         alpha_dual_ < 1);
}

void Ipm::makeStep() {
  stepSizes();

  // keep track of iterations with small stepsizes
  if (std::min(alpha_primal_, alpha_dual_) < 0.05)
    ++bad_iter_;
  else
    bad_iter_ = 0;

  // update iterate
  vectorAdd(it_->x, it_->delta.x, alpha_primal_);
  vectorAdd(it_->xl, it_->delta.xl, alpha_primal_);
  vectorAdd(it_->xu, it_->delta.xu, alpha_primal_);
  vectorAdd(it_->y, it_->delta.y, alpha_dual_);
  vectorAdd(it_->zl, it_->delta.zl, alpha_dual_);
  vectorAdd(it_->zu, it_->delta.zu, alpha_dual_);

  // compute new quantities
  it_->residual1234();
  it_->computeMu();
  it_->indicators();

  printOutput();
}

void Ipm::startingPoint() {
  std::vector<double>& x = it_->x;
  std::vector<double>& xl = it_->xl;
  std::vector<double>& xu = it_->xu;
  std::vector<double>& y = it_->y;
  std::vector<double>& zl = it_->zl;
  std::vector<double>& zu = it_->zu;

  // *********************************************************************
  // x starting point
  // *********************************************************************
  // compute feasible x
  for (Int i = 0; i < n_; ++i) {
    x[i] = 0.0;
    x[i] = std::max(x[i], model_.lb(i));
    x[i] = std::min(x[i], model_.ub(i));
  }

  const std::vector<double> temp_scaling(n_, 1.0);
  std::vector<double> temp_m(m_);

  if (options_.nla == kOptionNlaNormEq) {
    // use y to store b-A*x
    y = model_.b();
    model_.A().alphaProductPlusY(-1.0, x, y);

    // solve A*A^T * dx = b-A*x with factorisation and store the result in
    // temp_m

    // factorise A*A^T
    if (LS_->factorNE(model_.A(), temp_scaling)) goto failure;

    if (LS_->solveNE(y, temp_m)) goto failure;

  } else if (options_.nla == kOptionNlaAugmented) {
    // obtain solution of A*A^T * dx = b-A*x by solving
    // [ -I  A^T] [...] = [ -x]
    // [  A   0 ] [ dx] = [ b ]

    if (LS_->factorAS(model_.A(), temp_scaling)) goto failure;

    std::vector<double> rhs_x(n_);
    for (Int i = 0; i < n_; ++i) rhs_x[i] = -x[i];
    std::vector<double> lhs_x(n_);
    if (LS_->solveAS(rhs_x, model_.b(), lhs_x, temp_m)) goto failure;
  }

  // compute dx = A^T * (A*A^T)^{-1} * (b-A*x) and store the result in xl
  xl.assign(n_, 0.0);
  model_.A().alphaProductPlusY(1.0, temp_m, xl, true);

  // x += dx;
  vectorAdd(x, xl, 1.0);
  // *********************************************************************

  // *********************************************************************
  // xl, xu starting point
  // *********************************************************************
  // compute xl, xu that satisfy linear constraints
  {
    double violation{};
    for (Int i = 0; i < n_; ++i) {
      if (model_.hasLb(i)) {
        xl[i] = x[i] - model_.lb(i);
        violation = std::min(violation, xl[i]);
      } else {
        xl[i] = 0.0;
      }
      if (model_.hasUb(i)) {
        xu[i] = model_.ub(i) - x[i];
        violation = std::min(violation, xu[i]);
      } else {
        xu[i] = 0.0;
      }
    }

    // shift to be positive
    violation = 1.0 + std::max(0.0, -1.5 * violation);
    vectorAdd(xl, violation);
    vectorAdd(xu, violation);
  }
  // *********************************************************************

  // *********************************************************************
  // y starting point
  // *********************************************************************

  if (options_.nla == kOptionNlaNormEq) {
    // compute A*c
    std::fill(temp_m.begin(), temp_m.end(), 0.0);
    model_.A().alphaProductPlusY(1.0, model_.c(), temp_m);

    if (LS_->solveNE(temp_m, y)) goto failure;

  } else if (options_.nla == kOptionNlaAugmented) {
    // obtain solution of A*A^T * y = A*c by solving
    // [ -I  A^T] [...] = [ c ]
    // [  A   0 ] [ y ] = [ 0 ]

    std::vector<double> rhs_y(m_, 0.0);
    std::vector<double> lhs_x(n_);

    if (LS_->solveAS(model_.c(), rhs_y, lhs_x, y)) goto failure;
  }
  // *********************************************************************

  // *********************************************************************
  // zl, zu starting point
  // *********************************************************************
  // compute c - A^T * y and store in zl
  zl = model_.c();
  model_.A().alphaProductPlusY(-1.0, y, zl, true);

  // split result between zl and zu
  {
    double violation = 0.0;
    for (Int i = 0; i < n_; ++i) {
      double val = zl[i];
      zl[i] = 0.0;
      zu[i] = 0.0;

      if (model_.hasLb(i) && model_.hasUb(i)) {
        zl[i] = 0.5 * val;
        zu[i] = -0.5 * val;
      } else if (model_.hasLb(i)) {
        zl[i] = val;
      } else if (model_.hasUb(i)) {
        zu[i] = -val;
      }

      violation = std::min(violation, zl[i]);
      violation = std::min(violation, zu[i]);
    }

    // shift to be positive

    violation = 1.0 + std::max(0.0, -1.5 * violation);
    for (Int i = 0; i < n_; ++i) {
      if (model_.hasLb(i)) {
        zl[i] += violation;
      }
      if (model_.hasUb(i)) {
        zu[i] += violation;
      }
    }
  }
  // *********************************************************************

  // *********************************************************************
  // improve centrality
  // *********************************************************************
  {
    double xsum{1.0};
    double zsum{1.0};
    double mu{1.0};

    for (Int i = 0; i < n_; ++i) {
      if (model_.hasLb(i)) {
        xsum += xl[i];
        zsum += zl[i];
        mu += xl[i] * zl[i];
      }
      if (model_.hasUb(i)) {
        xsum += xu[i];
        zsum += zu[i];
        mu += xu[i] * zu[i];
      }
    }

    double dx = 0.5 * mu / zsum;
    double dz = 0.5 * mu / xsum;

    vectorAdd(xl, dx);
    vectorAdd(xu, dx);
    for (Int i = 0; i < n_; ++i) {
      if (model_.hasLb(i)) {
        zl[i] += dz;
      }
      if (model_.hasUb(i)) {
        zu[i] += dz;
      }
    }
  }
  // *********************************************************************

  return;

failure:
  std::cerr << "Error while computing starting point\n";
  info_.ipm_status = kIpmStatusError;
}

void Ipm::sigmaAffine() {
  sigma_ = kSigmaAffine;

  DataCollector::get()->setSigma(sigma_, true);
}

void Ipm::sigmaCorrectors() {
  if ((alpha_primal_ > 0.5 && alpha_dual_ > 0.5) || iter_ == 1) {
    sigma_ = 0.01;
  } else if (alpha_primal_ > 0.2 && alpha_dual_ > 0.2) {
    sigma_ = 0.1;
  } else if (alpha_primal_ > 0.1 && alpha_dual_ > 0.1) {
    sigma_ = 0.25;
  } else if (alpha_primal_ > 0.05 && alpha_dual_ > 0.05) {
    sigma_ = 0.5;
  } else {
    sigma_ = 0.9;
  }

  DataCollector::get()->setSigma(sigma_);
}

void Ipm::residualsMcc() {
  // compute right-hand side for multiple centrality correctors
  std::vector<double>& xl = it_->xl;
  std::vector<double>& xu = it_->xu;
  std::vector<double>& zl = it_->zl;
  std::vector<double>& zu = it_->zu;
  std::vector<double>& res5 = it_->res5;
  std::vector<double>& res6 = it_->res6;
  double& mu = it_->mu;

  // clear existing residuals
  it_->clearRes();

  // stepsizes of current direction
  double alpha_p, alpha_d;
  stepsToBoundary(alpha_p, alpha_d, it_->delta);

  // compute increased stepsizes
  alpha_p = std::max(1.0, alpha_p + kMccIncreaseAlpha);
  alpha_d = std::max(1.0, alpha_d + kMccIncreaseAlpha);

  // compute trial point
  std::vector<double> xlt = xl;
  std::vector<double> xut = xu;
  std::vector<double> zlt = zl;
  std::vector<double> zut = zu;
  vectorAdd(xlt, it_->delta.xl, alpha_p);
  vectorAdd(xut, it_->delta.xu, alpha_p);
  vectorAdd(zlt, it_->delta.zl, alpha_d);
  vectorAdd(zut, it_->delta.zu, alpha_d);

  // compute right-hand side for mcc
  for (Int i = 0; i < n_; ++i) {
    // res5
    if (model_.hasLb(i)) {
      double prod = xlt[i] * zlt[i];
      if (prod < sigma_ * mu * kGammaCorrector) {
        // prod is small, we add something positive to res5

        double temp = sigma_ * mu * kGammaCorrector - prod;
        res5[i] += temp;

      } else if (prod > sigma_ * mu / kGammaCorrector) {
        // prod is large, we may subtract something large from res5.
        // limit the amount to subtract to -sigma*mu/gamma

        double temp = sigma_ * mu / kGammaCorrector - prod;
        temp = std::max(temp, -sigma_ * mu / kGammaCorrector);
        res5[i] += temp;
      }
    } else {
      res5[i] = 0.0;
    }

    // res6
    if (model_.hasUb(i)) {
      double prod = xut[i] * zut[i];
      if (prod < sigma_ * mu * kGammaCorrector) {
        // prod is small, we add something positive to res6

        double temp = sigma_ * mu * kGammaCorrector - prod;
        res6[i] += temp;

      } else if (prod > sigma_ * mu / kGammaCorrector) {
        // prod is large, we may subtract something large from res6.
        // limit the amount to subtract to -sigma*mu/gamma

        double temp = sigma_ * mu / kGammaCorrector - prod;
        temp = std::max(temp, -sigma_ * mu / kGammaCorrector);
        res6[i] += temp;
      }
    } else {
      res6[i] = 0.0;
    }
  }
}

bool Ipm::centralityCorrectors() {
  // compute stepsizes of current direction
  double alpha_p_old, alpha_d_old;
  stepsToBoundary(alpha_p_old, alpha_d_old, it_->delta);

#ifdef PRINT_CORRECTORS
  printf("(%.2f,%.2f) -> ", alpha_p_old, alpha_d_old);
#endif

  Int cor;
  for (cor = 0; cor < info_.correctors; ++cor) {
    // compute rhs for corrector
    residualsMcc();

    // compute corrector
    NewtonDir corr(m_, n_);
    if (solveNewtonSystem(corr)) return true;
    if (recoverDirection(corr)) return true;

    double alpha_p, alpha_d;
    double wp = alpha_p_old * alpha_d_old;
    double wd = wp;
    bestWeight(it_->delta, corr, wp, wd, alpha_p, alpha_d);

#ifdef PRINT_CORRECTORS
    printf("(%.2f,%.2f) -> ", alpha_p, alpha_d);
#endif

    if (alpha_p < alpha_p_old + kMccIncreaseAlpha * kMccIncreaseMin &&
        alpha_d < alpha_d_old + kMccIncreaseAlpha * kMccIncreaseMin) {
      // reject corrector
#ifdef PRINT_CORRECTORS
      printf(" x");
#endif
      break;
    }

    if (alpha_p >= alpha_p_old + kMccIncreaseAlpha * kMccIncreaseMin) {
      // accept primal corrector
      vectorAdd(it_->delta.x, corr.x, wp);
      vectorAdd(it_->delta.xl, corr.xl, wp);
      vectorAdd(it_->delta.xu, corr.xu, wp);
      alpha_p_old = alpha_p;
    }
    if (alpha_d >= alpha_d_old + kMccIncreaseAlpha * kMccIncreaseMin) {
      // accept dual corrector
      vectorAdd(it_->delta.y, corr.y, wd);
      vectorAdd(it_->delta.zl, corr.zl, wd);
      vectorAdd(it_->delta.zu, corr.zu, wd);
      alpha_d_old = alpha_d;
    }

    if (alpha_p > 0.95 && alpha_d > 0.95) {
      // stepsizes are large enough, stop
      ++cor;
      break;
    }

    // else, keep computing correctors
  }
#ifdef PRINT_CORRECTORS
  printf("\n");
#endif

  DataCollector::get()->setCorrectors(cor);

  return false;
}

void Ipm::bestWeight(const NewtonDir& delta, const NewtonDir& corrector,
                     double& wp, double& wd, double& alpha_p,
                     double& alpha_d) const {
  // Find the best primal and dual weights for the corrector in the interval
  // [alpha_p_old * alpha_d_old, 1].
  // Upon return, wp and wd are the optimal weights, alpha_p and alpha_d are the
  // corresponding stepsizes.

  // keep track of best stepsizes
  alpha_p = 0.0;
  alpha_d = 0.0;

  // initial weight
  double w = wp;

  // divide interval into 9 points
  const double step = (1.0 - w) / 8;

  // for each weight, compute stepsizes and save best ones
  for (; w <= 1.0; w += step) {
    double ap, ad;
    stepsToBoundary(ap, ad, delta, &corrector, w);
    if (ap > alpha_p) {
      alpha_p = ap;
      wp = w;
    }
    if (ad > alpha_d) {
      alpha_d = ad;
      wd = w;
    }

    if (step == 0.0) break;
  }
}

bool Ipm::checkIterate() {
  // Check that iterate is not NaN or Inf
  if (it_->isNan()) {
    printf("\nIterate is nan\n");
    info_.ipm_status = kIpmStatusError;
    return true;
  } else if (it_->isInf()) {
    printf("\nIterate is inf\n");
    info_.ipm_status = kIpmStatusError;
    return true;
  }

  // check that no component is negative
  for (Int i = 0; i < n_; ++i) {
    if ((model_.hasLb(i) && it_->xl[i] < 0) ||
        (model_.hasLb(i) && it_->zl[i] < 0) ||
        (model_.hasUb(i) && it_->xu[i] < 0) ||
        (model_.hasUb(i) && it_->zu[i] < 0)) {
      printf("\nIterate has negative component\n");
      return true;
    }
  }

  return false;
}

bool Ipm::checkBadIter() {
  // If too many bad iterations, stop
  if (bad_iter_ >= kMaxBadIter) {
    printf("=== No progress\n");
    info_.ipm_status = kIpmStatusNoProgress;
    return true;
  }
  return false;
}

bool Ipm::checkTermination() {
  bool feasible = it_->pinf < options_.feasibility_tol &&
                  it_->dinf < options_.feasibility_tol;
  bool optimal = it_->pdgap < options_.optimality_tol;

  bool terminate = false;

  if (feasible && optimal) {
    if (info_.ipm_status != kIpmStatusPDFeas)
      printf("=== Primal-dual feasible point found\n");

    info_.ipm_status = kIpmStatusPDFeas;

    if (options_.crossover) {
      bool ready_for_crossover =
          it_->infeasAfterDropping() < options_.crossover_tol;
      if (ready_for_crossover) {
        printf("=== Ready for crossover\n");
        terminate = true;
      }
    } else {
      terminate = true;
    }
  }
  return terminate;
}

void Ipm::backwardError(const NewtonDir& delta) const {
  std::vector<double>& x = it_->x;
  std::vector<double>& xl = it_->xl;
  std::vector<double>& xu = it_->xu;
  std::vector<double>& y = it_->y;
  std::vector<double>& zl = it_->zl;
  std::vector<double>& zu = it_->zu;
  std::vector<double>& res1 = it_->res1;
  std::vector<double>& res2 = it_->res2;
  std::vector<double>& res3 = it_->res3;
  std::vector<double>& res4 = it_->res4;
  std::vector<double>& res5 = it_->res5;
  std::vector<double>& res6 = it_->res6;

  // ===================================================================================
  // Normwise backward error
  // ===================================================================================

  // residuals of the six blocks of equations
  // res1 - A * dx
  std::vector<double> r1 = res1;
  model_.A().alphaProductPlusY(-1.0, delta.x, r1);

  // res2 - dx + dxl
  std::vector<double> r2(n_);
  for (Int i = 0; i < n_; ++i) r2[i] = res2[i] - delta.x[i] + delta.xl[i];

  // res3 - dx - dxu
  std::vector<double> r3(n_);
  for (Int i = 0; i < n_; ++i) r3[i] = res3[i] - delta.x[i] - delta.xu[i];

  // res4 - A^T * dy - dzl + dzu
  std::vector<double> r4(n_);
  for (Int i = 0; i < n_; ++i) r4[i] = res4[i] - delta.zl[i] + delta.zu[i];
  model_.A().alphaProductPlusY(-1.0, delta.y, r4, true);

  // res5 - Zl * Dxl - Xl * Dzl
  std::vector<double> r5(n_);
  for (Int i = 0; i < n_; ++i) {
    if (model_.hasLb(i))
      r5[i] = res5[i] - zl[i] * delta.xl[i] - xl[i] * delta.zl[i];
  }

  // res6 - Zu * Dxu - Xu * Dzu
  std::vector<double> r6(n_);
  for (Int i = 0; i < n_; ++i) {
    if (model_.hasUb(i))
      r6[i] = res6[i] - zu[i] * delta.xu[i] - xu[i] * delta.zu[i];
  }

  // ...and their infinity norm
  double inf_norm_r{};
  inf_norm_r = std::max(inf_norm_r, infNorm(r1));
  inf_norm_r = std::max(inf_norm_r, infNorm(r2));
  inf_norm_r = std::max(inf_norm_r, infNorm(r3));
  inf_norm_r = std::max(inf_norm_r, infNorm(r4));
  inf_norm_r = std::max(inf_norm_r, infNorm(r5));
  inf_norm_r = std::max(inf_norm_r, infNorm(r6));

  // infinity norm of solution
  double inf_norm_delta{};
  inf_norm_delta = std::max(inf_norm_delta, infNorm(delta.x));
  inf_norm_delta = std::max(inf_norm_delta, infNorm(delta.xl));
  inf_norm_delta = std::max(inf_norm_delta, infNorm(delta.xu));
  inf_norm_delta = std::max(inf_norm_delta, infNorm(delta.y));
  inf_norm_delta = std::max(inf_norm_delta, infNorm(delta.zl));
  inf_norm_delta = std::max(inf_norm_delta, infNorm(delta.zu));

  // infinity norm of rhs
  double inf_norm_res{};
  inf_norm_res = std::max(inf_norm_res, infNorm(res1));
  inf_norm_res = std::max(inf_norm_res, infNorm(res2));
  inf_norm_res = std::max(inf_norm_res, infNorm(res3));
  inf_norm_res = std::max(inf_norm_res, infNorm(res4));
  inf_norm_res = std::max(inf_norm_res, infNorm(res5));
  inf_norm_res = std::max(inf_norm_res, infNorm(res6));

  // infinity norm of big 6x6 matrix:
  // max( ||A||_inf, 2, 2+||A||_1, max_j(zl_j+xl_j), max_j(zu_j+xu_j) )

  std::vector<double> norm_cols_A(n_);
  std::vector<double> norm_rows_A(m_);
  for (Int col = 0; col < n_; ++col) {
    for (Int el = model_.A().start_[col]; el < model_.A().start_[col + 1];
         ++el) {
      Int row = model_.A().index_[el];
      double val = model_.A().value_[el];
      norm_cols_A[col] += std::abs(val);
      norm_rows_A[row] += std::abs(val);
    }
  }
  double one_norm_A = *std::max_element(norm_cols_A.begin(), norm_cols_A.end());
  double inf_norm_A = *std::max_element(norm_rows_A.begin(), norm_rows_A.end());

  double inf_norm_matrix = inf_norm_A;
  inf_norm_matrix = std::max(inf_norm_matrix, one_norm_A + 2);
  for (Int i = 0; i < n_; ++i) {
    if (model_.hasLb(i))
      inf_norm_matrix = std::max(inf_norm_matrix, zl[i] + xl[i]);
    if (model_.hasUb(i))
      inf_norm_matrix = std::max(inf_norm_matrix, zu[i] + xu[i]);
  }

  // compute normwise backward error:
  // ||residual|| / ( ||matrix|| * ||solution|| + ||rhs|| )
  double nw_back_err =
      inf_norm_r / (inf_norm_matrix * inf_norm_delta + inf_norm_res);

  // ===================================================================================
  // Componentwise backward error
  // ===================================================================================

  // Compute |A| * |dx| and |A^T| * |dy|
  std::vector<double> abs_prod_A(m_);
  std::vector<double> abs_prod_At(n_);
  for (Int col = 0; col < n_; ++col) {
    for (Int el = model_.A().start_[col]; el < model_.A().start_[col + 1];
         ++el) {
      Int row = model_.A().index_[el];
      double val = model_.A().value_[el];
      abs_prod_A[row] += std::abs(val) * std::abs(delta.x[col]);
      abs_prod_At[col] += std::abs(val) * std::abs(delta.y[row]);
    }
  }

  // componentwise backward error:
  // max |residual_i| / (|matrix| * |solution| + |rhs|)_i
  double cw_back_err{};

  // first block
  for (Int i = 0; i < m_; ++i) {
    double denom = abs_prod_A[i] + std::abs(res1[i]);
    double num = std::abs(r1[i]);
    if (denom == 0.0) {
      if (num != 0.0) cw_back_err = std::numeric_limits<double>::max();
    } else
      cw_back_err = std::max(cw_back_err, num / denom);
  }
  // second and third block
  for (Int i = 0; i < n_; ++i) {
    if (model_.hasLb(i)) {
      double denom =
          std::abs(delta.x[i]) + std::abs(delta.xl[i]) + std::abs(res2[i]);
      double num = std::abs(r2[i]);
      if (denom == 0.0) {
        if (num != 0.0) cw_back_err = std::numeric_limits<double>::max();
      } else
        cw_back_err = std::max(cw_back_err, num / denom);
    }
    if (model_.hasUb(i)) {
      double denom =
          std::abs(delta.x[i]) + std::abs(delta.xu[i]) + std::abs(res3[i]);
      double num = std::abs(r3[i]);
      if (denom == 0.0) {
        if (num != 0.0) cw_back_err = std::numeric_limits<double>::max();
      } else
        cw_back_err = std::max(cw_back_err, num / denom);
    }
  }
  // fourth block
  for (Int i = 0; i < n_; ++i) {
    double denom = abs_prod_At[i] + std::abs(res4[i]);
    if (model_.hasLb(i)) denom += std::abs(delta.zl[i]);
    if (model_.hasUb(i)) denom += std::abs(delta.zu[i]);
    double num = std::abs(r4[i]);
    if (denom == 0.0) {
      if (num != 0.0) cw_back_err = std::numeric_limits<double>::max();
    } else
      cw_back_err = std::max(cw_back_err, num / denom);
  }
  // fifth and sixth block
  for (Int i = 0; i < n_; ++i) {
    if (model_.hasLb(i)) {
      double denom = zl[i] * std::abs(delta.xl[i]) +
                     xl[i] * std::abs(delta.zl[i]) + std::abs(res5[i]);
      double num = std::abs(r5[i]);
      if (denom == 0.0) {
        if (num != 0.0) cw_back_err = std::numeric_limits<double>::max();
      } else
        cw_back_err = std::max(cw_back_err, num / denom);
    }
    if (model_.hasUb(i)) {
      double denom = zu[i] * std::abs(delta.xu[i]) +
                     xu[i] * std::abs(delta.zu[i]) + std::abs(res6[i]);
      double num = std::abs(r6[i]);
      if (denom == 0.0) {
        if (num != 0.0) cw_back_err = std::numeric_limits<double>::max();
      } else
        cw_back_err = std::max(cw_back_err, num / denom);
    }
  }

  DataCollector::get()->setBackError(nw_back_err, cw_back_err);
}

void Ipm::printHeader() const {
  if (iter_ % 20 == 0) {
    printf(
        " iter      primal obj        dual obj        pinf      dinf "
        "       mu      alpha p/d    p/d gap    time\n");
  }
}

void Ipm::printOutput() const {
  printHeader();

  printf(
      "%5d %16.8e %16.8e %10.2e %10.2e %10.2e %6.2f %5.2f %9.2e "
      "%7.1f\n",
      iter_, it_->pobj, it_->dobj, it_->pinf, it_->dinf, it_->mu, alpha_primal_,
      alpha_dual_, it_->pdgap, clock_.stop());
}

void Ipm::printInfo() const {
  printf("\n");
  printf("Problem %s\n", model_.name().c_str());
  printf("%.2e rows, %.2e cols, %.2e nnz\n", (double)m_, (double)n_,
         (double)model_.A().numNz());

  if (options_.parallel == kOptionParallelOff)
    printf("Running on 1 thread\n");
  else
    printf("Running on %d threads\n", highs::parallel::num_threads());

#ifdef DEBUG
  printf("Running in debug mode\n");
#endif

  printf("\n");

  // print range of coefficients
  model_.checkCoefficients();
}

const IpmInfo& Ipm::getInfo() const { return info_; }
void Ipm::getSolution(std::vector<double>& x, std::vector<double>& xl,
                      std::vector<double>& xu, std::vector<double>& slack,
                      std::vector<double>& y, std::vector<double>& zl,
                      std::vector<double>& zu) const {
  // prepare and return solution with internal format

  if (!info_.ipx_used) {
    it_->extract(x, xl, xu, slack, y, zl, zu);
    model_.unscale(x, xl, xu, slack, y, zl, zu);
    model_.postprocess(slack, y);
  } else {
    ipx_lps_.GetInteriorSolution(x.data(), xl.data(), xu.data(), slack.data(),
                                 y.data(), zl.data(), zu.data());
  }
}

void Ipm::getSolution(std::vector<double>& x, std::vector<double>& slack,
                      std::vector<double>& y, std::vector<double>& z) const {
  // prepare and return solution with format for crossover

  if (!info_.ipx_used) {
    it_->extract(x, slack, y, z);
    model_.unscale(x, slack, y, z);
    model_.postprocess(slack, y);
  } else {
    ipx_lps_.GetBasicSolution(x.data(), slack.data(), y.data(), z.data(),
                              nullptr, nullptr);
  }
}

void Ipm::maxCorrectors() {
  if (kMaxCorrectors > 0) {
    // Compute estimate of effort to factorise and solve

    // Effort to factorise depends on the number of flops
    double fact_effort = LS_->flops() + 100 * LS_->spops();

    // Effort to solve depends on the number of nonzeros of L multiplied by 2,
    // because there are two sweeps through L (forward and backward).
    double solv_effort = 2.0 * LS_->nz();

    // The factorise phase uses BLAS-3 and can be parallelised, the solve phase
    // uses BLAS-2 and cannot be parallelised. To account for this, the
    // factorisation effort is multiplied by a coefficient < 1, estimated
    // empirically.
    double alpha = 1.0 / 112.0;

    double ratio = alpha * fact_effort / solv_effort;

    // At each ipm iteration, there are up to (1+k) directions computed, where k
    // is the number of correctors. Each direction requires up (1+f) solves,
    // where f is the number of refinement steps. So, up to (1+k)(1+f) solves
    // are performed per iteration. However, not all refinement steps are used
    // all the time, so use f/2.
    // Therefore, we want (1+k)(1+f/2) < ratio.

    double thresh = ratio / (1.0 + kMaxRefinementIter / 2.0) - 1;

    info_.correctors = std::floor(thresh);
    info_.correctors = std::max(info_.correctors, (Int)1);
    info_.correctors = std::min(info_.correctors, kMaxCorrectors);

  } else {
    info_.correctors = -kMaxCorrectors;
  }
}

}  // namespace highspm