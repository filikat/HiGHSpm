#include "Ipm.h"

#include <cassert>
#include <cmath>
#include <iostream>

#include "parallel/HighsParallel.h"

void Ipm::load(const int num_var, const int num_con, const double* obj,
               const double* rhs, const double* lower, const double* upper,
               const int* A_ptr, const int* A_rows, const double* A_vals,
               const int* constraints, const std::string& pb_name,
               const Options& options) {
  if (!obj || !rhs || !lower || !upper || !A_ptr || !A_rows || !A_vals ||
      !constraints)
    return;

  model_.init(num_var, num_con, obj, rhs, lower, upper, A_ptr, A_rows, A_vals,
              constraints, pb_name);

  model_.scale();
  model_.reformulate();

  m_ = model_.num_con_;
  n_ = model_.num_var_;

  options_ = options;
}

Output Ipm::solve() {
  if (!model_.ready_) return Output{};

  printInfo();

  // ------------------------------------------
  // ---- INITIALIZE --------------------------
  // ------------------------------------------

  // start timer
  clock_.start();

  // initialize iterate and residuals
  it_ = Iterate(m_, n_);
  res_ = Residuals(m_, n_);

  DataCollector::start();

  // initialize linear solver
  LS_.reset(new FactorHiGHSSolver(options_));
  if (LS_->setup(model_.A_, options_)) return Output{};
  LS_->clear();

  // initialize starting point, residuals and mu
  computeStartingPoint();
  computeResiduals1234();
  computeMu();
  computeIndicators();

  // ------------------------------------------
  // ---- MAIN LOOP ---------------------------
  // ------------------------------------------

  while (iter_ < kMaxIterations) {
    if (checkIterate()) break;
    if (checkBadIter()) break;
    if (checkTermination()) break;

    ++iter_;

    // Clear Newton direction
    delta_ = NewtonDir(m_, n_);

    // Clear any existing data in the linear solver
    LS_->clear();

    computeScaling();

    // ===== PREDICTOR =====
    sigma_ = kSigmaAffine;
    computeResiduals56();
    if (solveNewtonSystem(delta_)) break;
    if (recoverDirection(delta_)) break;

    // ===== CORRECTORS =====
    computeSigma();
    if (centralityCorrectors()) break;

    // ===== STEP =====
    makeStep();
    computeResiduals1234();
    computeMu();
    computeIndicators();
    collectData();
    printOutput();
  }

  LS_->finalise();
  model_.unscale(it_);

  // output struct
  Output out{};
  out.it = std::move(it_);
  out.iterations = iter_;
  out.primal_infeas = primal_infeas_;
  out.dual_infeas = dual_infeas_;
  out.mu = mu_;
  out.status = ipm_status_;

  DataCollector::get()->printIter();
  DataCollector::destruct();

  return out;
}

void Ipm::computeMu() {
  mu_ = 0.0;
  int number_finite_bounds{};
  for (int i = 0; i < n_; ++i) {
    if (model_.hasLb(i)) {
      mu_ += it_.xl[i] * it_.zl[i];
      ++number_finite_bounds;
    }
    if (model_.hasUb(i)) {
      mu_ += it_.xu[i] * it_.zu[i];
      ++number_finite_bounds;
    }
  }
  mu_ /= number_finite_bounds;
}

void Ipm::computeResiduals1234() {
  // res1
  res_.res1 = model_.b_;
  model_.A_.alphaProductPlusY(-1.0, it_.x, res_.res1);

  // res2
  for (int i = 0; i < n_; ++i) {
    if (model_.hasLb(i)) {
      res_.res2[i] = model_.lower_[i] - it_.x[i] + it_.xl[i];
    } else {
      res_.res2[i] = 0.0;
    }
  }

  // res3
  for (int i = 0; i < n_; ++i) {
    if (model_.hasUb(i)) {
      res_.res3[i] = model_.upper_[i] - it_.x[i] - it_.xu[i];
    } else {
      res_.res3[i] = 0.0;
    }
  }

  // res4
  res_.res4 = model_.c_;
  model_.A_.alphaProductPlusY(-1.0, it_.y, res_.res4, true);
  for (int i = 0; i < n_; ++i) {
    if (model_.hasLb(i)) {
      res_.res4[i] -= it_.zl[i];
    }
    if (model_.hasUb(i)) {
      res_.res4[i] += it_.zu[i];
    }
  }

  // Check for NaN or Inf
  assert(!res_.isNaN());
  assert(!res_.isInf());
}

void Ipm::computeResiduals56() {
  // For predictor, sigma_ should be zero (or small).
  // For corrector, delta_ should be the affine scaling direction and sigma_
  // should be computed with computeSigma().

  for (int i = 0; i < n_; ++i) {
    // res5
    if (model_.hasLb(i)) {
      res_.res5[i] = sigma_ * mu_ - it_.xl[i] * it_.zl[i];
    } else {
      res_.res5[i] = 0.0;
    }

    // res6
    if (model_.hasUb(i)) {
      res_.res6[i] = sigma_ * mu_ - it_.xu[i] * it_.zu[i];
    } else {
      res_.res6[i] = 0.0;
    }
  }
}

std::vector<double> Ipm::computeResiduals7() {
  std::vector<double> res7;

  res7 = res_.res4;
  for (int i = 0; i < n_; ++i) {
    if (model_.hasLb(i)) {
      res7[i] -= ((res_.res5[i] + it_.zl[i] * res_.res2[i]) / it_.xl[i]);
    }
    if (model_.hasUb(i)) {
      res7[i] += ((res_.res6[i] - it_.zu[i] * res_.res3[i]) / it_.xu[i]);
    }
  }

  return res7;
}

std::vector<double> Ipm::computeResiduals8(const std::vector<double>& res7) {
  std::vector<double> res8;

  res8 = res_.res1;

  std::vector<double> temp(res7);

  // temp = (Theta^-1+Rp)^-1 * res7
  for (int i = 0; i < n_; ++i) {
    temp[i] /= scaling_[i] + kPrimalStaticRegularization;
  }

  // res8 += A * temp
  model_.A_.alphaProductPlusY(1.0, temp, res8);

  return res8;
}

void Ipm::computeScaling() {
  scaling_.assign(n_, 0.0);

  for (int i = 0; i < n_; ++i) {
    if (model_.hasLb(i)) {
      scaling_[i] += it_.zl[i] / it_.xl[i];
    }
    if (model_.hasUb(i)) {
      scaling_[i] += it_.zu[i] / it_.xu[i];
    }
    // slow down the growth of theta
    if (scaling_[i] < 1e-12) scaling_[i] = sqrt(1e-12 * scaling_[i]);
  }
}

bool Ipm::solveNewtonSystem(NewtonDir& delta) {
  std::vector<double> res7{computeResiduals7()};

  // NORMAL EQUATIONS
  if (options_.nla == kOptionNlaNormEq) {
    std::vector<double> res8{computeResiduals8(res7)};

    // factorise normal equations, if not yet done
    if (!LS_->valid_ && LS_->factorNE(model_.A_, scaling_)) goto failure;

    // solve with normal equations
    if (LS_->solveNE(res8, delta.y)) goto failure;

    // Compute delta.x
    // Deltax = A^T * Deltay - res7;
    delta.x = res7;
    model_.A_.alphaProductPlusY(-1.0, delta.y, delta.x, true);
    vectorScale(delta.x, -1.0);

    // Deltax = (Theta^-1+Rp)^-1 * Deltax
    for (int i = 0; i < n_; ++i)
      delta.x[i] /= scaling_[i] + kPrimalStaticRegularization;

  }

  // AUGMENTED SYSTEM
  else {
    // factorise augmented system, if not yet done
    if (!LS_->valid_ && LS_->factorAS(model_.A_, scaling_)) goto failure;

    // solve with augmented system
    if (LS_->solveAS(res7, res_.res1, delta.x, delta.y)) goto failure;
  }

  // iterative refinement
  LS_->refine(model_.A_, scaling_, res7, res_.res1, delta.x, delta.y);

  return false;

// Failure occured in factorisation or solve
failure:
  std::cerr << "Error while solving Newton system\n";
  ipm_status_ = "Error";
  return true;
}

bool Ipm::recoverDirection(NewtonDir& delta) {
  // Recover components xl, xu, zl, zu of partial direction delta.

  for (int i = 0; i < n_; ++i) {
    if (model_.hasLb(i) || model_.hasUb(i)) {
      delta.xl[i] = delta.x[i] - res_.res2[i];
      delta.zl[i] = (res_.res5[i] - it_.zl[i] * delta.xl[i]) / it_.xl[i];
    } else {
      delta.xl[i] = 0.0;
      delta.zl[i] = 0.0;
    }
  }
  for (int i = 0; i < n_; ++i) {
    if (model_.hasLb(i) || model_.hasUb(i)) {
      delta.xu[i] = res_.res3[i] - delta.x[i];
      delta.zu[i] = (res_.res6[i] - it_.zu[i] * delta.xu[i]) / it_.xu[i];
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
      if (std::isfinite(it_.xl[i]) && std::isfinite(it_.xu[i])) {
        if (it_.zl[i] * it_.xu[i] >= it_.zu[i] * it_.xl[i])
          delta.zl[i] = res_.res4[i] + delta.zu[i] - Atdy[i];
        else
          delta.zu[i] = -res_.res4[i] + delta.zl[i] + Atdy[i];
      } else if (std::isfinite(it_.xl[i])) {
        delta.zl[i] = res_.res4[i] + delta.zu[i] - Atdy[i];
      } else {
        delta.zu[i] = -res_.res4[i] + delta.zl[i] + Atdy[i];
      }
    }
  }

  // Check for NaN of Inf
  if (delta.isNaN()) {
    std::cerr << "Direction is nan\n";
    ipm_status_ = "Error";
    return true;
  } else if (delta.isInf()) {
    std::cerr << "Direciton is inf\n";
    ipm_status_ = "Error";
    return true;
  }
  return false;
}

void Ipm::computeStepSizes(const NewtonDir& delta, const NewtonDir* corr,
                           double& alpha_primal, double& alpha_dual) const {
  // Compute primal and dual stepsizes, based on direction delta.
  // If corr is valid, then the stepsizes are computed based on direction
  // delta+corr.

  alpha_primal = 1.0;
  double p_limit_x = 0.0;
  double p_limit_dx = 0.0;
  for (int i = 0; i < n_; ++i) {
    double add = (corr ? corr->xl[i] : 0.0);
    if ((delta.xl[i] + add) < 0 && model_.hasLb(i)) {
      double val = -it_.xl[i] / (delta.xl[i] + add);
      if (val < alpha_primal) {
        alpha_primal = val;
        p_limit_x = it_.xl[i];
        p_limit_dx = -(delta.xl[i] + add);
      }
    }
    add = (corr ? corr->xu[i] : 0.0);
    if ((delta.xu[i] + add) < 0 && model_.hasUb(i)) {
      double val = -it_.xu[i] / (delta.xu[i] + add);
      if (val < alpha_primal) {
        alpha_primal = val;
        p_limit_x = it_.xu[i];
        p_limit_dx = -(delta.xu[i] + add);
      }
    }
  }
  alpha_primal *= kInteriorScaling;

  alpha_dual = 1.0;
  double d_limit_z = 0.0;
  double d_limit_dz = 0.0;
  for (int i = 0; i < n_; ++i) {
    double add = (corr ? corr->zl[i] : 0.0);
    if ((delta.zl[i] + add) < 0 && model_.hasLb(i)) {
      double val = -it_.zl[i] / (delta.zl[i] + add);
      if (val < alpha_dual) {
        alpha_dual = val;
        d_limit_z = it_.zl[i];
        d_limit_dz = -(delta.zl[i] + add);
      }
    }
    add = (corr ? corr->zu[i] : 0.0);
    if ((delta.zu[i] + add) < 0 && model_.hasUb(i)) {
      double val = -it_.zu[i] / (delta.zu[i] + add);
      if (val < alpha_dual) {
        alpha_dual = val;
        d_limit_z = it_.zu[i];
        d_limit_dz = -(delta.zu[i] + add);
      }
    }
  }
  alpha_dual *= kInteriorScaling;

  DataCollector::get()->back().p_limit_x = p_limit_x;
  DataCollector::get()->back().p_limit_dx = p_limit_dx;
  DataCollector::get()->back().d_limit_z = d_limit_z;
  DataCollector::get()->back().d_limit_dz = d_limit_dz;

  assert(alpha_primal > 0 && alpha_primal < 1 && alpha_dual > 0 &&
         alpha_dual < 1);
}

void Ipm::makeStep() {
  computeStepSizes(delta_, nullptr, alpha_primal_, alpha_dual_);

  if (std::min(alpha_primal_, alpha_dual_) < 0.05)
    ++bad_iter_;
  else
    bad_iter_ = 0;

  vectorAdd(it_.x, delta_.x, alpha_primal_);
  vectorAdd(it_.xl, delta_.xl, alpha_primal_);
  vectorAdd(it_.xu, delta_.xu, alpha_primal_);
  vectorAdd(it_.y, delta_.y, alpha_dual_);
  vectorAdd(it_.zl, delta_.zl, alpha_dual_);
  vectorAdd(it_.zu, delta_.zu, alpha_dual_);
}

void Ipm::computeStartingPoint() {
  // use conjugate gradient for starting point
  CgSolver CG;

  // *********************************************************************
  // x starting point
  // *********************************************************************
  // compute feasible x
  for (int i = 0; i < n_; ++i) {
    it_.x[i] = 0.0;
    it_.x[i] = std::max(it_.x[i], model_.lower_[i]);
    it_.x[i] = std::min(it_.x[i], model_.upper_[i]);
  }

  const std::vector<double> temp_scaling(n_, 1.0);
  std::vector<double> temp_m(m_, 1.0);

  // use y to store b-A*x
  it_.y = model_.b_;
  model_.A_.alphaProductPlusY(-1.0, it_.x, it_.y);

  // solve A*A^T * dx = b-A*x with factorization and store the result in
  // temp_m

  // factorize A*A^T
  CG.factorNE(model_.A_, temp_scaling);

  int itercg1 = CG.solveNE(it_.y, temp_m);

  // compute dx = A^T * (A*A^T)^{-1} * (b-A*x) and store the result in xl
  it_.xl.assign(n_, 0.0);
  model_.A_.alphaProductPlusY(1.0, temp_m, it_.xl, true);

  // x += dx;
  vectorAdd(it_.x, it_.xl, 1.0);
  // *********************************************************************

  // *********************************************************************
  // xl, xu starting point
  // *********************************************************************
  // compute xl, xu that satisfy linear constraints
  {
    double violation{};
    for (int i = 0; i < n_; ++i) {
      if (model_.hasLb(i)) {
        it_.xl[i] = it_.x[i] - model_.lower_[i];
        violation = std::min(violation, it_.xl[i]);
      } else {
        it_.xl[i] = 0.0;
      }
      if (model_.hasUb(i)) {
        it_.xu[i] = model_.upper_[i] - it_.x[i];
        violation = std::min(violation, it_.xu[i]);
      } else {
        it_.xu[i] = 0.0;
      }
    }

    // shift to be positive
    violation = 1.0 + std::max(0.0, -1.5 * violation);
    vectorAdd(it_.xl, violation);
    vectorAdd(it_.xu, violation);
  }
  // *********************************************************************

  // *********************************************************************
  // y starting point
  // *********************************************************************

  // compute A*c
  temp_m.assign(m_, 0.0);
  model_.A_.alphaProductPlusY(1.0, model_.c_, temp_m);

  int itercg2 = CG.solveNE(temp_m, it_.y);

  // *********************************************************************

  // *********************************************************************
  // zl, zu starting point
  // *********************************************************************
  // compute c - A^T * y and store in zl
  it_.zl = model_.c_;
  model_.A_.alphaProductPlusY(-1.0, it_.y, it_.zl, true);

  // split result between zl and zu
  {
    double violation = 0.0;
    for (int i = 0; i < n_; ++i) {
      double val = it_.zl[i];
      it_.zl[i] = 0.0;
      it_.zu[i] = 0.0;

      if (model_.hasLb(i) && model_.hasUb(i)) {
        it_.zl[i] = 0.5 * val;
        it_.zu[i] = -0.5 * val;
      } else if (model_.hasLb(i)) {
        it_.zl[i] = val;
      } else if (model_.hasUb(i)) {
        it_.zu[i] = -val;
      }

      violation = std::min(violation, it_.zl[i]);
      violation = std::min(violation, it_.zu[i]);
    }

    // shift to be positive

    violation = 1.0 + std::max(0.0, -1.5 * violation);
    for (int i = 0; i < n_; ++i) {
      if (model_.hasLb(i)) {
        it_.zl[i] += violation;
      }
      if (model_.hasUb(i)) {
        it_.zu[i] += violation;
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
        xsum += it_.xl[i];
        zsum += it_.zl[i];
        mu += it_.xl[i] * it_.zl[i];
      }
      if (model_.hasUb(i)) {
        xsum += it_.xu[i];
        zsum += it_.zu[i];
        mu += it_.xu[i] * it_.zu[i];
      }
    }

    double dx = 0.5 * mu / zsum;
    double dz = 0.5 * mu / xsum;

    vectorAdd(it_.xl, dx);
    vectorAdd(it_.xu, dx);
    for (int i = 0; i < n_; ++i) {
      if (model_.hasLb(i)) {
        it_.zl[i] += dz;
      }
      if (model_.hasUb(i)) {
        it_.zu[i] += dz;
      }
    }
  }
  // *********************************************************************

  printf("Starting point required %d + %d CG iterations\n\n", itercg1, itercg2);
}

void Ipm::computeSigma() {
  /*if (min_prod_ < kSmallProduct || max_prod_ > kLargeProduct) {
    // bad complementarity products, perform centring
    sigma_ = 0.9;
  } else*/
  // good complementarity products, decide based on previous iteration
  if ((alpha_primal_ > 0.5 && alpha_dual_ > 0.5) || iter_ == 1) {
    sigma_ = 0.01;
  } else if (alpha_primal_ > 0.1 && alpha_dual_ > 0.1) {
    sigma_ = 0.1;
  } else if (alpha_primal_ > 0.05 && alpha_dual_ > 0.05) {
    sigma_ = 0.25;
  } else if (alpha_primal_ > 0.02 && alpha_dual_ > 0.02) {
    sigma_ = 0.5;
  } else {
    sigma_ = 0.9;
  }

  DataCollector::get()->back().sigma = sigma_;
}

void Ipm::computeResidualsMcc() {
  // compute right-hand side for multiple centrality correctors

  // clear existing residuals
  res_ = Residuals(m_, n_);

  // stepsizes of current direction
  double alpha_p, alpha_d;
  computeStepSizes(delta_, nullptr, alpha_p, alpha_d);

  // compute increased stepsizes
  alpha_p = std::max(1.0, alpha_p + kMccIncreaseAlpha);
  alpha_d = std::max(1.0, alpha_d + kMccIncreaseAlpha);

  // compute trial point
  std::vector<double> xlt = it_.xl;
  std::vector<double> xut = it_.xu;
  std::vector<double> zlt = it_.zl;
  std::vector<double> zut = it_.zu;
  vectorAdd(xlt, delta_.xl, alpha_p);
  vectorAdd(xut, delta_.xu, alpha_p);
  vectorAdd(zlt, delta_.zl, alpha_d);
  vectorAdd(zut, delta_.zu, alpha_d);

  double max_v1{};
  double max_v2{};

  // compute right-hand side for mcc
  for (int i = 0; i < n_; ++i) {
    // res5
    if (model_.hasLb(i)) {
      double prod = xlt[i] * zlt[i];
      if (prod < sigma_ * mu_ * kGammaCorrector) {
        // prod is small, we add something positive to res5

        double temp = sigma_ * mu_ * kGammaCorrector - prod;
        res_.res5[i] += temp;
        max_v1 = std::max(max_v1, temp);

      } else if (prod > sigma_ * mu_ / kGammaCorrector) {
        // prod is large, we may subtract something large from res5.
        // limit the amount to subtract to -sigma*mu/gamma

        double temp = sigma_ * mu_ / kGammaCorrector - prod;
        temp = std::max(temp, -sigma_ * mu_ / kGammaCorrector);
        res_.res5[i] += temp;
        max_v2 = std::max(max_v2, std::abs(temp));
      }
    } else {
      res_.res5[i] = 0.0;
    }

    // res6
    if (model_.hasUb(i)) {
      double prod = xut[i] * zut[i];
      if (prod < sigma_ * mu_ * kGammaCorrector) {
        // prod is small, we add something positive to res6

        double temp = sigma_ * mu_ * kGammaCorrector - prod;
        res_.res6[i] += temp;
        max_v1 = std::max(max_v1, temp);

      } else if (prod > sigma_ * mu_ / kGammaCorrector) {
        // prod is large, we may subtract something large from res6.
        // limit the amount to subtract to -sigma*mu/gamma

        double temp = sigma_ * mu_ / kGammaCorrector - prod;
        temp = std::max(temp, -sigma_ * mu_ / kGammaCorrector);
        res_.res6[i] += temp;
        max_v2 = std::max(max_v2, std::abs(temp));
      }
    } else {
      res_.res6[i] = 0.0;
    }
  }

  // printf("%e %e\n", max_v1, max_v2);
}

bool Ipm::centralityCorrectors() {
  // compute stepsizes of current direction
  double alpha_p_old, alpha_d_old;
  computeStepSizes(delta_, nullptr, alpha_p_old, alpha_d_old);

#ifdef PRINT_CORRECTORS
  printf("(%.2f,%.2f) -> ", alpha_p_old, alpha_d_old);
#endif

  int cor;
  for (cor = 0; cor < kMaxCorrectors; ++cor) {
    // compute rhs for corrector
    computeResidualsMcc();

    // compute corrector
    NewtonDir corr(m_, n_);
    if (solveNewtonSystem(corr)) return true;
    if (recoverDirection(corr)) return true;

    // stepsizes of new corrected direction
    double alpha_p, alpha_d;
    computeStepSizes(delta_, &corr, alpha_p, alpha_d);

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
      vectorAdd(delta_.x, corr.x);
      vectorAdd(delta_.xl, corr.xl);
      vectorAdd(delta_.xu, corr.xu);
      alpha_p_old = alpha_p;
    }
    if (alpha_d >= alpha_d_old + kMccIncreaseAlpha * kMccIncreaseMin) {
      // accept dual corrector
      vectorAdd(delta_.y, corr.y);
      vectorAdd(delta_.zl, corr.zl);
      vectorAdd(delta_.zu, corr.zu);
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

void Ipm::computeIndicators() {
  primal_infeas_ = infNorm(res_.res1);
  primal_infeas_ = std::max(primal_infeas_, infNorm(res_.res2));
  primal_infeas_ = std::max(primal_infeas_, infNorm(res_.res3));
  primal_infeas_ /= (1 + model_.normRhs());

  dual_infeas_ = infNorm(res_.res4) / (1 + model_.normObj());

  // compute unscaled primal and dual objectives
  primal_obj_ = dotProd(it_.x, model_.c_);

  if (model_.colexp_.size() > 0)
    primal_obj_ = std::ldexp(primal_obj_, -model_.cexp_ - model_.bexp_);

  dual_obj_ = dotProd(it_.y, model_.b_);
  for (int i = 0; i < n_; ++i) {
    if (model_.hasLb(i)) {
      dual_obj_ += model_.lower_[i] * it_.zl[i];
    }
    if (model_.hasUb(i)) {
      dual_obj_ -= model_.upper_[i] * it_.zu[i];
    }
  }

  if (model_.colexp_.size() > 0)
    dual_obj_ = std::ldexp(dual_obj_, -model_.cexp_ - model_.bexp_);

  pd_gap_ = std::abs(primal_obj_ - dual_obj_) /
            (1 + 0.5 * std::abs(primal_obj_ + dual_obj_));

  computeProducts();
}

void Ipm::computeProducts() {
  if (iter_ == 0) return;

  // compute min and max entry in Theta
  double& min_theta = DataCollector::get()->back().min_theta;
  double& max_theta = DataCollector::get()->back().max_theta;
  min_theta = kInf;
  max_theta = 0.0;
  for (int i = 0; i < n_; ++i) {
    if (scaling_[i] != 0.0) {
      min_theta = std::min(min_theta, 1.0 / scaling_[i]);
      max_theta = std::max(max_theta, 1.0 / scaling_[i]);
    }
  }

  // compute min and max complementarity product
  // (x_l)_j * (z_l)_j / mu or (x_u)_j * (z_u)_j / mu

  min_prod_ = std::numeric_limits<double>::max();
  max_prod_ = 0.0;
  int& num_small = DataCollector::get()->back().num_small_prod;
  int& num_large = DataCollector::get()->back().num_large_prod;

  for (int i = 0; i < n_; ++i) {
    if (model_.hasLb(i)) {
      double prod = it_.xl[i] * it_.zl[i] / mu_;
      min_prod_ = std::min(min_prod_, prod);
      max_prod_ = std::max(max_prod_, prod);
      if (prod < kSmallProduct) {
        ++num_small;
        /*printf("Small product %e %e\n", it_.xl[i], it_.zl[i]);
        double ratio = it_.xl[i] / it_.zl[i];
        if (ratio < 1.0 / kThreshProduct) {
          // xl[i] is too small, perturb it
          it_.xl[i] = kSmallProduct * mu_ / it_.zl[i];
          printf("\txl set to %e\n", it_.xl[i]);
        } else if (ratio > kThreshProduct) {
          // zl[i] is too small, perturb it
          it_.zl[i] = kSmallProduct * mu_ / it_.xl[i];
          printf("\tzl set to %e\n", it_.zl[i]);
        }*/
      }
      if (prod > kLargeProduct) {
        ++num_large;
        /*printf("Large product %e %e\n", it_.xl[i], it_.zl[i]);
        double ratio = it_.xl[i] / it_.zl[i];
        if (ratio < 1.0 / kThreshProduct) {
          // zl[i] is too large, perturb it
          it_.zl[i] = kLargeProduct * mu_ / it_.xl[i];
          printf("\tzl set to %e\n", it_.zl[i]);
        } else if (ratio > kThreshProduct) {
          // xl[i] is too large, perturb it
          it_.xl[i] = kLargeProduct * mu_ / it_.zl[i];
          printf("\txl set to %e\n", it_.xl[i]);
        }*/
      }
    }
    if (model_.hasUb(i)) {
      double prod = it_.xu[i] * it_.zu[i] / mu_;
      min_prod_ = std::min(min_prod_, prod);
      max_prod_ = std::max(max_prod_, prod);
      if (prod < kSmallProduct) {
        ++num_small;
        /*printf("Small product %e %e\n", it_.xu[i], it_.zu[i]);
        double ratio = it_.xu[i] / it_.zu[i];
        if (ratio < 1.0 / kThreshProduct) {
          // xu[i] is too small, perturb it
          it_.xu[i] = kSmallProduct * mu_ / it_.zu[i];
          printf("\txu set to %e\n", it_.xu[i]);
        } else if (ratio > kThreshProduct) {
          // zu[i] is too small, perturb it
          it_.zu[i] = kSmallProduct * mu_ / it_.xu[i];
          printf("\tzu set to %e\n", it_.zu[i]);
        }*/
      }
      if (prod > kLargeProduct) {
        ++num_large;
        /*printf("Large product %e %e\n", it_.xu[i], it_.zu[i]);
        double ratio = it_.xu[i] / it_.zu[i];
        if (ratio < 1.0 / kThreshProduct) {
          // zu[i] is too large, perturb it
          it_.zu[i] = kLargeProduct * mu_ / it_.xu[i];
          printf("\tzu set to %e\n", it_.zu[i]);
        } else if (ratio > kThreshProduct) {
          // xu[i] is too large, perturb it
          it_.xu[i] = kLargeProduct * mu_ / it_.zu[i];
          printf("\txu set to %e\n", it_.xu[i]);
        }*/
      }
    }
  }
  DataCollector::get()->back().min_prod = min_prod_;
  DataCollector::get()->back().max_prod = max_prod_;
}

bool Ipm::checkIterate() {
  // Check that iterate is not NaN or Inf
  if (it_.isNaN()) {
    std::cerr << "iterate is nan\n";
    ipm_status_ = "Error";
    return true;
  } else if (it_.isInf()) {
    std::cerr << "iterate is inf\n";
    ipm_status_ = "Error";
    return true;
  }
  return false;
}

bool Ipm::checkBadIter() {
  // If too many bad iterations, stop
  if (bad_iter_ >= kMaxBadIter) {
    printf("\n Failure: no progress\n\n");
    ipm_status_ = "No progress";
    return true;
  }
  return false;
}

bool Ipm::checkTermination() {
  if (pd_gap_ < kIpmTolerance &&         // primal-dual gap is small
      primal_infeas_ < kIpmTolerance &&  // primal feasibility
      dual_infeas_ < kIpmTolerance) {    // dual feasibility
    printf("\n===== Optimal solution found =====\n");

    ipm_status_ = "Optimal";
    return true;
  }
  return false;
}

void Ipm::printHeader() const {
  if (iter_ % 20 == 1) {
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
      iter_, primal_obj_, dual_obj_, primal_infeas_, dual_infeas_, mu_,
      alpha_primal_, alpha_dual_, pd_gap_, clock_.stop());
}

void Ipm::printInfo() const {
  printf("\n");
  printf("Problem %s\n", model_.pb_name_.c_str());
  printf("%.2e rows, %.2e cols, %.2e nnz\n", (double)m_, (double)n_,
         (double)model_.A_.numNz());
  printf("Using %s\n", options_.nla == kOptionNlaAugmented
                           ? "augmented systems"
                           : "normal equations");

#ifdef PARALLEL_TREE
  printf("Running on %d threads\n", highs::parallel::num_threads());
#else
  printf("Running on 1 thread\n");
#endif

  printf("\n");

  // print range of coefficients
  model_.checkCoefficients();
}

void Ipm::collectData() const {
  DataCollector::get()->back().p_obj = primal_obj_;
  DataCollector::get()->back().d_obj = dual_obj_;
  DataCollector::get()->back().p_inf = primal_infeas_;
  DataCollector::get()->back().d_inf = dual_infeas_;
  DataCollector::get()->back().mu = mu_;
  DataCollector::get()->back().pd_gap = pd_gap_;
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
      minxl = std::min(minxl, std::abs(it_.xl[i]));
      maxxl = std::max(maxxl, std::abs(it_.xl[i]));
      minzl = std::min(minzl, std::abs(it_.zl[i]));
      maxzl = std::max(maxzl, std::abs(it_.zl[i]));
      mindxl = std::min(mindxl, std::abs(delta_.xl[i]));
      maxdxl = std::max(maxdxl, std::abs(delta_.xl[i]));
      mindzl = std::min(mindzl, std::abs(delta_.zl[i]));
      maxdzl = std::max(maxdzl, std::abs(delta_.zl[i]));
    }
    if (model_.hasUb(i)) {
      minxu = std::min(minxu, std::abs(it_.xu[i]));
      maxxu = std::max(maxxu, std::abs(it_.xu[i]));
      minzu = std::min(minzu, std::abs(it_.zu[i]));
      maxzu = std::max(maxzu, std::abs(it_.zu[i]));
      mindxu = std::min(mindxu, std::abs(delta_.xu[i]));
      maxdxu = std::max(maxdxu, std::abs(delta_.xu[i]));
      mindzu = std::min(mindzu, std::abs(delta_.zu[i]));
      maxdzu = std::max(maxdzu, std::abs(delta_.zu[i]));
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