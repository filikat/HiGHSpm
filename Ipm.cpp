#include "Ipm.h"

#include <cassert>
#include <cmath>
#include <iostream>

#include "parallel/HighsParallel.h"

void Ipm::load(const int num_var, const int num_con, const double* obj,
               const double* rhs, const double* lower, const double* upper,
               const int* A_ptr, const int* A_rows, const double* A_vals,
               const char* constraints, const std::string& pb_name,
               const Options& options) {
  if (!obj || !rhs || !lower || !upper || !A_ptr || !A_rows || !A_vals ||
      !constraints)
    return;

  model_.init(num_var, num_con, obj, rhs, lower, upper, A_ptr, A_rows, A_vals,
              constraints, pb_name);

  model_.scale();
  model_.reformulate();

  m_ = model_.m_;
  n_ = model_.n_;

  options_ = options;
}

IpmStatus Ipm::solve() {
  if (!model_.ready_) return kIpmStatusError;

  printInfo();

  // ------------------------------------------
  // ---- INITIALIZE --------------------------
  // ------------------------------------------

  // start timer
  clock_.start();

  DataCollector::start();

  // initialize iterate object
  it_.reset(new IpmIterate(model_));

  // initialize linear solver
  LS_.reset(new FactorHiGHSSolver(options_));
  if (LS_->setup(model_.A_, options_)) return kIpmStatusError;
  LS_->clear();

  // initialize starting point, residuals and mu
  startingPoint();
  it_->residual1234();
  it_->mu();
  it_->indicators();
  printOutput();

  // ------------------------------------------
  // ---- MAIN LOOP ---------------------------
  // ------------------------------------------

  while (iter_ < kMaxIterations) {
    if (checkIterate()) break;
    if (checkBadIter()) break;
    if (checkTermination()) break;

    ++iter_;

    // Clear Newton direction
    it_->clearDir();

    // Clear any existing data in the linear solver
    LS_->clear();

    // compute theta inverse
    it_->scaling();

    // ===== PREDICTOR =====
    sigmaAffine();
    it_->residual56(sigma_);
    if (solveNewtonSystem(it_->delta_)) break;
    if (recoverDirection(it_->delta_)) break;

    // ===== CORRECTORS =====
    sigmaCorrectors();
    if (centralityCorrectors()) break;

    // ===== STEP =====
    makeStep();
    it_->residual1234();
    it_->mu();
    it_->indicators();
    collectData();
    printOutput();
  }

  LS_->finalise();

  // solution for the user
  it_->prepareForUser(x_user, xl_user, xu_user, slack_user, y_user, zl_user,
                      zu_user);

  DataCollector::get()->printIter();
  DataCollector::destruct();

  return ipm_status_;
}

bool Ipm::solveNewtonSystem(NewtonDir& delta) {
  std::vector<double>& theta_inv = it_->scaling_;

  std::vector<double> res7 = it_->residual7();

  // NORMAL EQUATIONS
  if (options_.nla == kOptionNlaNormEq) {
    std::vector<double> res8 = it_->residual8(res7);

    // factorise normal equations, if not yet done
    if (!LS_->valid_ && LS_->factorNE(model_.A_, theta_inv)) goto failure;

    // solve with normal equations
    if (LS_->solveNE(res8, delta.y)) goto failure;

    // Compute delta.x
    // Deltax = A^T * Deltay - res7;
    delta.x = res7;
    model_.A_.alphaProductPlusY(-1.0, delta.y, delta.x, true);
    vectorScale(delta.x, -1.0);

    // Deltax = (Theta^-1+Rp)^-1 * Deltax
    for (int i = 0; i < n_; ++i)
      delta.x[i] /= theta_inv[i] + kPrimalStaticRegularization;

  }

  // AUGMENTED SYSTEM
  else {
    // factorise augmented system, if not yet done
    if (!LS_->valid_ && LS_->factorAS(model_.A_, theta_inv)) goto failure;

    // solve with augmented system
    if (LS_->solveAS(res7, it_->res1_, delta.x, delta.y)) goto failure;
  }

  return false;

// Failure occured in factorisation or solve
failure:
  std::cerr << "Error while solving Newton system\n";
  ipm_status_ = kIpmStatusError;
  return true;
}

bool Ipm::recoverDirection(NewtonDir& delta) {
  // Recover components xl, xu, zl, zu of partial direction delta.
  std::vector<double>& xl = it_->xl_;
  std::vector<double>& xu = it_->xu_;
  std::vector<double>& zl = it_->zl_;
  std::vector<double>& zu = it_->zu_;
  std::vector<double>& res2 = it_->res2_;
  std::vector<double>& res3 = it_->res3_;
  std::vector<double>& res4 = it_->res4_;
  std::vector<double>& res5 = it_->res5_;
  std::vector<double>& res6 = it_->res6_;

  for (int i = 0; i < n_; ++i) {
    if (model_.hasLb(i) || model_.hasUb(i)) {
      delta.xl[i] = delta.x[i] - res2[i];
      delta.zl[i] = (res5[i] - zl[i] * delta.xl[i]) / xl[i];
    } else {
      delta.xl[i] = 0.0;
      delta.zl[i] = 0.0;
    }
  }
  for (int i = 0; i < n_; ++i) {
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
  model_.A_.alphaProductPlusY(1.0, delta.y, Atdy, true);
  for (int i = 0; i < n_; ++i) {
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
    ipm_status_ = kIpmStatusError;
    return true;
  } else if (it_->isDirInf()) {
    std::cerr << "Direciton is inf\n";
    ipm_status_ = kIpmStatusError;
    return true;
  }
  return false;
}

double Ipm::stepToBoundary(const std::vector<double>& x,
                           const std::vector<double>& dx,
                           const std::vector<double>* cor, double weight,
                           bool lo, int* block) const {
  // Compute the largest alpha s.t. x + alpha * dx >= 0.
  // If cor is valid, consider x + alpha * (dx + w * cor) instead.
  // Use lo=1 for xl and zl, lo=0 for xu and zu.
  // Return the blocking index in block.

  const double damp = 1.0 - std::numeric_limits<double>::epsilon();

  double alpha = 1.0;
  int bl = -1;

  for (int i = 0; i < x.size(); ++i) {
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

  int block;
  double axl = stepToBoundary(it_->xl_, delta.xl, cor ? &(cor->xl) : nullptr,
                              weight, true);
  double axu = stepToBoundary(it_->xu_, delta.xu, cor ? &(cor->xu) : nullptr,
                              weight, false);
  double azl = stepToBoundary(it_->zl_, delta.zl, cor ? &(cor->zl) : nullptr,
                              weight, true);
  double azu = stepToBoundary(it_->zu_, delta.zu, cor ? &(cor->zu) : nullptr,
                              weight, false);

  alpha_primal = std::min(axl, axu);
  alpha_primal = std::min(alpha_primal, 1.0);
  alpha_dual = std::min(azl, azu);
  alpha_dual = std::min(alpha_dual, 1.0);
}

void Ipm::stepSizes() {
  // Compute primal and dual stepsizes.
  std::vector<double>& xl = it_->xl_;
  std::vector<double>& xu = it_->xu_;
  std::vector<double>& zl = it_->zl_;
  std::vector<double>& zu = it_->zu_;
  std::vector<double>& dxl = it_->dxl_;
  std::vector<double>& dxu = it_->dxu_;
  std::vector<double>& dzl = it_->dzl_;
  std::vector<double>& dzu = it_->dzu_;

  // parameters for Mehrotra heuristic
  const double gamma_f = 0.9;
  const double gamma_a = 1.0 / (1.0 - gamma_f);

  // compute stepsizes and blocking components
  int block_xl, block_xu, block_zl, block_zu;
  double alpha_xl = stepToBoundary(xl, dxl, nullptr, 0, true, &block_xl);
  double alpha_xu = stepToBoundary(xu, dxu, nullptr, 0, false, &block_xu);
  double alpha_zl = stepToBoundary(zl, dzl, nullptr, 0, true, &block_zl);
  double alpha_zu = stepToBoundary(zu, dzu, nullptr, 0, false, &block_zu);

  double max_p = std::min(alpha_xl, alpha_xu);
  double max_d = std::min(alpha_zl, alpha_zu);

  // compute mu with current stepsizes
  double mu_full = 0.0;
  int num_finite = 0;
  for (int i = 0; i < n_; ++i) {
    if (model_.hasLb(i)) {
      mu_full +=
          (xl[i] + max_p * dxl[i]) * (zl[i] + max_d * dzl[i]);
      ++num_finite;
    }
    if (model_.hasUb(i)) {
      mu_full +=
          (xu[i] + max_p * dxu[i]) * (zu[i] + max_d * dzu[i]);
      ++num_finite;
    }
  }
  mu_full /= num_finite;
  mu_full /= gamma_a;

  // compute new stepsizes based on Mehrotra heuristic

  // primal
  double alpha_p = 1.0;
  int block_p = -1;
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
  int block_d = -1;
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

  if (std::min(alpha_primal_, alpha_dual_) < 0.05)
    ++bad_iter_;
  else
    bad_iter_ = 0;

  vectorAdd(it_->x_, it_->dx_, alpha_primal_);
  vectorAdd(it_->xl_, it_->dxl_, alpha_primal_);
  vectorAdd(it_->xu_, it_->dxu_, alpha_primal_);
  vectorAdd(it_->y_, it_->dy_, alpha_dual_);
  vectorAdd(it_->zl_, it_->dzl_, alpha_dual_);
  vectorAdd(it_->zu_, it_->dzu_, alpha_dual_);
}

void Ipm::startingPoint() {
  std::vector<double>& x = it_->x_;
  std::vector<double>& xl = it_->xl_;
  std::vector<double>& xu = it_->xu_;
  std::vector<double>& y = it_->y_;
  std::vector<double>& zl = it_->zl_;
  std::vector<double>& zu = it_->zu_;

  // *********************************************************************
  // x starting point
  // *********************************************************************
  // compute feasible x
  for (int i = 0; i < n_; ++i) {
    x[i] = 0.0;
    x[i] = std::max(x[i], model_.lower_[i]);
    x[i] = std::min(x[i], model_.upper_[i]);
  }

  const std::vector<double> temp_scaling(n_, 1.0);
  std::vector<double> temp_m(m_);

  if (options_.nla == kOptionNlaNormEq) {
    // use y to store b-A*x
    y = model_.b_;
    model_.A_.alphaProductPlusY(-1.0, x, y);

    // solve A*A^T * dx = b-A*x with factorization and store the result in
    // temp_m

    // factorize A*A^T
    if (LS_->factorNE(model_.A_, temp_scaling)) goto failure;

    if (LS_->solveNE(y, temp_m)) goto failure;

  } else if (options_.nla == kOptionNlaAugmented) {
    // obtain solution of A*A^T * dx = b-A*x by solving
    // [ -I  A^T] [...] = [ -x]
    // [  A   0 ] [ dx] = [ b ]

    if (LS_->factorAS(model_.A_, temp_scaling)) goto failure;

    std::vector<double> rhs_x(n_);
    for (int i = 0; i < n_; ++i) rhs_x[i] = -x[i];
    std::vector<double> lhs_x(n_);
    if (LS_->solveAS(rhs_x, model_.b_, lhs_x, temp_m)) goto failure;
  }

  // compute dx = A^T * (A*A^T)^{-1} * (b-A*x) and store the result in xl
  xl.assign(n_, 0.0);
  model_.A_.alphaProductPlusY(1.0, temp_m, xl, true);

  // x += dx;
  vectorAdd(x, xl, 1.0);
  // *********************************************************************

  // *********************************************************************
  // xl, xu starting point
  // *********************************************************************
  // compute xl, xu that satisfy linear constraints
  {
    double violation{};
    for (int i = 0; i < n_; ++i) {
      if (model_.hasLb(i)) {
        xl[i] = x[i] - model_.lower_[i];
        violation = std::min(violation, xl[i]);
      } else {
        xl[i] = 0.0;
      }
      if (model_.hasUb(i)) {
        xu[i] = model_.upper_[i] - x[i];
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
    model_.A_.alphaProductPlusY(1.0, model_.c_, temp_m);

    if (LS_->solveNE(temp_m, y)) goto failure;

  } else if (options_.nla == kOptionNlaAugmented) {
    // obtain solution of A*A^T * y = A*c by solving
    // [ -I  A^T] [...] = [ c ]
    // [  A   0 ] [ y ] = [ 0 ]

    std::vector<double> rhs_y(m_, 0.0);
    std::vector<double> lhs_x(n_);

    if (LS_->solveAS(model_.c_, rhs_y, lhs_x, y)) goto failure;
  }
  // *********************************************************************

  // *********************************************************************
  // zl, zu starting point
  // *********************************************************************
  // compute c - A^T * y and store in zl
  zl = model_.c_;
  model_.A_.alphaProductPlusY(-1.0, y, zl, true);

  // split result between zl and zu
  {
    double violation = 0.0;
    for (int i = 0; i < n_; ++i) {
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
    for (int i = 0; i < n_; ++i) {
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

    for (int i = 0; i < n_; ++i) {
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
    for (int i = 0; i < n_; ++i) {
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
  ipm_status_ = kIpmStatusError;
}

void Ipm::sigmaAffine() {
  sigma_ = kSigmaAffine;

  DataCollector::get()->back().sigma_aff = sigma_;
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

  DataCollector::get()->back().sigma = sigma_;
}

void Ipm::residualsMcc() {
  // compute right-hand side for multiple centrality correctors
  std::vector<double>& xl = it_->xl_;
  std::vector<double>& xu = it_->xu_;
  std::vector<double>& zl = it_->zl_;
  std::vector<double>& zu = it_->zu_;
  std::vector<double>& res5 = it_->res5_;
  std::vector<double>& res6 = it_->res6_;
  double& mu = it_->mu_;

  // clear existing residuals
  it_->clearRes();

  // stepsizes of current direction
  double alpha_p, alpha_d;
  stepsToBoundary(alpha_p, alpha_d, it_->delta_);

  // compute increased stepsizes
  alpha_p = std::max(1.0, alpha_p + kMccIncreaseAlpha);
  alpha_d = std::max(1.0, alpha_d + kMccIncreaseAlpha);

  // compute trial point
  std::vector<double> xlt = xl;
  std::vector<double> xut = xu;
  std::vector<double> zlt = zl;
  std::vector<double> zut = zu;
  vectorAdd(xlt, it_->dxl_, alpha_p);
  vectorAdd(xut, it_->dxu_, alpha_p);
  vectorAdd(zlt, it_->dzl_, alpha_d);
  vectorAdd(zut, it_->dzu_, alpha_d);

  // compute right-hand side for mcc
  for (int i = 0; i < n_; ++i) {
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
  stepsToBoundary(alpha_p_old, alpha_d_old, it_->delta_);

#ifdef PRINT_CORRECTORS
  printf("(%.2f,%.2f) -> ", alpha_p_old, alpha_d_old);
#endif

  int cor;
  for (cor = 0; cor < kMaxCorrectors; ++cor) {
    // compute rhs for corrector
    residualsMcc();

    // compute corrector
    NewtonDir corr(m_, n_);
    if (solveNewtonSystem(corr)) return true;
    if (recoverDirection(corr)) return true;

    double alpha_p, alpha_d;
    double wp = alpha_p_old * alpha_d_old;
    double wd = wp;
    bestWeight(it_->delta_, corr, wp, wd, alpha_p, alpha_d);

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
      vectorAdd(it_->dx_, corr.x, wp);
      vectorAdd(it_->dxl_, corr.xl, wp);
      vectorAdd(it_->dxu_, corr.xu, wp);
      alpha_p_old = alpha_p;
    }
    if (alpha_d >= alpha_d_old + kMccIncreaseAlpha * kMccIncreaseMin) {
      // accept dual corrector
      vectorAdd(it_->dy_, corr.y, wd);
      vectorAdd(it_->dzl_, corr.zl, wd);
      vectorAdd(it_->dzu_, corr.zu, wd);
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

  DataCollector::get()->back().correctors = cor;

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
    ipm_status_ = kIpmStatusError;
    return true;
  } else if (it_->isInf()) {
    printf("\nIterate is inf\n");
    ipm_status_ = kIpmStatusError;
    return true;
  }

  // check that no component is negative
  for (int i = 0; i < n_; ++i) {
    if ((model_.hasLb(i) && it_->xl_[i] < 0) ||
        (model_.hasLb(i) && it_->zl_[i] < 0) ||
        (model_.hasUb(i) && it_->xu_[i] < 0) ||
        (model_.hasUb(i) && it_->zu_[i] < 0)) {
      printf("\nIterative has negative component\n");
      return true;
    }
  }

  return false;
}

bool Ipm::checkBadIter() {
  // If too many bad iterations, stop
  if (bad_iter_ >= kMaxBadIter) {
    printf("\n Failure: no progress\n\n");
    ipm_status_ = kIpmStatusNoProgress;
    return true;
  }
  return false;
}

bool Ipm::checkTermination() {
  if (it_->pdgap_ < kIpmTolerance &&  // primal-dual gap is small
      it_->pinf_ < kIpmTolerance &&   // primal feasibility
      it_->dinf_ < kIpmTolerance) {   // dual feasibility
    printf("\n===== Optimal solution found =====\n\n");

    ipm_status_ = kIpmStatusOptimal;
    return true;
  }
  return false;
}

void Ipm::backwardError(const NewtonDir& delta) const {
  std::vector<double>& x = it_->x_;
  std::vector<double>& xl = it_->xl_;
  std::vector<double>& xu = it_->xu_;
  std::vector<double>& y = it_->y_;
  std::vector<double>& zl = it_->zl_;
  std::vector<double>& zu = it_->zu_;
  std::vector<double>& res1 = it_->res1_;
  std::vector<double>& res2 = it_->res2_;
  std::vector<double>& res3 = it_->res3_;
  std::vector<double>& res4 = it_->res4_;
  std::vector<double>& res5 = it_->res5_;
  std::vector<double>& res6 = it_->res6_;

  // ===================================================================================
  // Normwise backward error
  // ===================================================================================

  // residuals of the six blocks of equations
  // res1 - A * dx
  std::vector<double> r1 = res1;
  model_.A_.alphaProductPlusY(-1.0, delta.x, r1);

  // res2 - dx + dxl
  std::vector<double> r2(n_);
  for (int i = 0; i < n_; ++i) r2[i] = res2[i] - delta.x[i] + delta.xl[i];

  // res3 - dx - dxu
  std::vector<double> r3(n_);
  for (int i = 0; i < n_; ++i) r3[i] = res3[i] - delta.x[i] - delta.xu[i];

  // res4 - A^T * dy - dzl + dzu
  std::vector<double> r4(n_);
  for (int i = 0; i < n_; ++i) r4[i] = res4[i] - delta.zl[i] + delta.zu[i];
  model_.A_.alphaProductPlusY(-1.0, delta.y, r4, true);

  // res5 - Zl * Dxl - Xl * Dzl
  std::vector<double> r5(n_);
  for (int i = 0; i < n_; ++i) {
    if (model_.hasLb(i))
      r5[i] = res5[i] - zl[i] * delta.xl[i] - xl[i] * delta.zl[i];
  }

  // res6 - Zu * Dxu - Xu * Dzu
  std::vector<double> r6(n_);
  for (int i = 0; i < n_; ++i) {
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
  for (int col = 0; col < n_; ++col) {
    for (int el = model_.A_.start_[col]; el < model_.A_.start_[col + 1]; ++el) {
      int row = model_.A_.index_[el];
      double val = model_.A_.value_[el];
      norm_cols_A[col] += std::abs(val);
      norm_rows_A[row] += std::abs(val);
    }
  }
  double one_norm_A = *std::max_element(norm_cols_A.begin(), norm_cols_A.end());
  double inf_norm_A = *std::max_element(norm_rows_A.begin(), norm_rows_A.end());

  double inf_norm_matrix = inf_norm_A;
  inf_norm_matrix = std::max(inf_norm_matrix, one_norm_A + 2);
  for (int i = 0; i < n_; ++i) {
    if (model_.hasLb(i))
      inf_norm_matrix = std::max(inf_norm_matrix, zl[i] + xl[i]);
    if (model_.hasUb(i))
      inf_norm_matrix = std::max(inf_norm_matrix, zu[i] + xu[i]);
  }

  // compute normwise backward error:
  // ||residual|| / ( ||matrix|| * ||solution|| + ||rhs|| )
  double nw_back_err =
      inf_norm_r / (inf_norm_matrix * inf_norm_delta + inf_norm_res);

  DataCollector::get()->back().nw_back_err =
      std::max(DataCollector::get()->back().nw_back_err, nw_back_err);

  // ===================================================================================
  // Componentwise backward error
  // ===================================================================================

  // Compute |A| * |dx| and |A^T| * |dy|
  std::vector<double> abs_prod_A(m_);
  std::vector<double> abs_prod_At(n_);
  for (int col = 0; col < n_; ++col) {
    for (int el = model_.A_.start_[col]; el < model_.A_.start_[col + 1]; ++el) {
      int row = model_.A_.index_[el];
      double val = model_.A_.value_[el];
      abs_prod_A[row] += std::abs(val) * std::abs(delta.x[col]);
      abs_prod_At[col] += std::abs(val) * std::abs(delta.y[row]);
    }
  }

  // componentwise backward error:
  // max |residual_i| / (|matrix| * |solution| + |rhs|)_i
  double cw_back_err{};

  // first block
  for (int i = 0; i < m_; ++i) {
    double denom = abs_prod_A[i] + std::abs(res1[i]);
    double num = std::abs(r1[i]);
    if (denom == 0.0) {
      if (num != 0.0) cw_back_err = std::numeric_limits<double>::max();
    } else
      cw_back_err = std::max(cw_back_err, num / denom);
  }
  // second and third block
  for (int i = 0; i < n_; ++i) {
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
  for (int i = 0; i < n_; ++i) {
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
  for (int i = 0; i < n_; ++i) {
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

  DataCollector::get()->back().cw_back_err =
      std::max(DataCollector::get()->back().cw_back_err, cw_back_err);
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
      iter_, it_->pobj_, it_->dobj_, it_->pinf_, it_->dinf_, it_->mu_,
      alpha_primal_, alpha_dual_, it_->pdgap_, clock_.stop());
}

void Ipm::printInfo() const {
  printf("\n");
  printf("Problem %s\n", model_.pb_name_.c_str());
  printf("%.2e rows, %.2e cols, %.2e nnz\n", (double)m_, (double)n_,
         (double)model_.A_.numNz());
  printf("Using %s\n", options_.nla == kOptionNlaAugmented
                           ? "augmented systems"
                           : "normal equations");

#if (defined(PARALLEL_TREE) || defined(PARALLEL_NODE))
  printf("Running on %d threads\n", highs::parallel::num_threads());
#else
  printf("Running on 1 thread\n");
#endif

  printf("\n");

  // print range of coefficients
  model_.checkCoefficients();
}

void Ipm::collectData() const {
  DataCollector::get()->back().p_obj = it_->pobj_;
  DataCollector::get()->back().d_obj = it_->dobj_;
  DataCollector::get()->back().p_inf = it_->pinf_;
  DataCollector::get()->back().d_inf = it_->dinf_;
  DataCollector::get()->back().mu = it_->mu_;
  DataCollector::get()->back().pd_gap = it_->pdgap_;
  DataCollector::get()->back().p_alpha = alpha_primal_;
  DataCollector::get()->back().d_alpha = alpha_dual_;

  double& minxl = DataCollector::get()->back().min_xl;
  double& minxu = DataCollector::get()->back().min_xu;
  double& minzl = DataCollector::get()->back().min_zl;
  double& minzu = DataCollector::get()->back().min_zu;
  double& maxxl = DataCollector::get()->back().max_xl;
  double& maxxu = DataCollector::get()->back().max_xu;
  double& maxzl = DataCollector::get()->back().max_zl;
  double& maxzu = DataCollector::get()->back().max_zu;

  double& mindxl = DataCollector::get()->back().min_dxl;
  double& mindxu = DataCollector::get()->back().min_dxu;
  double& mindzl = DataCollector::get()->back().min_dzl;
  double& mindzu = DataCollector::get()->back().min_dzu;
  double& maxdxl = DataCollector::get()->back().max_dxl;
  double& maxdxu = DataCollector::get()->back().max_dxu;
  double& maxdzl = DataCollector::get()->back().max_dzl;
  double& maxdzu = DataCollector::get()->back().max_dzu;

  for (int i = 0; i < n_; ++i) {
    if (model_.hasLb(i)) {
      minxl = std::min(minxl, std::abs(it_->xl_[i]));
      maxxl = std::max(maxxl, std::abs(it_->xl_[i]));
      minzl = std::min(minzl, std::abs(it_->zl_[i]));
      maxzl = std::max(maxzl, std::abs(it_->zl_[i]));
      mindxl = std::min(mindxl, std::abs(it_->dxl_[i]));
      maxdxl = std::max(maxdxl, std::abs(it_->dxl_[i]));
      mindzl = std::min(mindzl, std::abs(it_->dzl_[i]));
      maxdzl = std::max(maxdzl, std::abs(it_->dzl_[i]));
    }
    if (model_.hasUb(i)) {
      minxu = std::min(minxu, std::abs(it_->xu_[i]));
      maxxu = std::max(maxxu, std::abs(it_->xu_[i]));
      minzu = std::min(minzu, std::abs(it_->zu_[i]));
      maxzu = std::max(maxzu, std::abs(it_->zu_[i]));
      mindxu = std::min(mindxu, std::abs(it_->dxu_[i]));
      maxdxu = std::max(maxdxu, std::abs(it_->dxu_[i]));
      mindzu = std::min(mindzu, std::abs(it_->dzu_[i]));
      maxdzu = std::max(maxdzu, std::abs(it_->dzu_[i]));
    }
  }

  if (minxl == std::numeric_limits<double>::max()) minxl = 0.0;
  if (minxu == std::numeric_limits<double>::max()) minxu = 0.0;
  if (minzl == std::numeric_limits<double>::max()) minzl = 0.0;
  if (minzu == std::numeric_limits<double>::max()) minzu = 0.0;
  if (mindxl == std::numeric_limits<double>::max()) mindxl = 0.0;
  if (mindxu == std::numeric_limits<double>::max()) mindxu = 0.0;
  if (mindzl == std::numeric_limits<double>::max()) mindzl = 0.0;
  if (mindzu == std::numeric_limits<double>::max()) mindzu = 0.0;
}

int Ipm::getIter() const { return iter_; }
void Ipm::getSolution(std::vector<double>& x, std::vector<double>& xl,
                      std::vector<double>& xu, std::vector<double>& slack,
                      std::vector<double>& y, std::vector<double>& zl,
                      std::vector<double>& zu) const {
  x = x_user;
  xl = xl_user;
  xu = xu_user;
  slack = slack_user;
  y = y_user;
  zl = zl_user;
  zu = zu_user;
}