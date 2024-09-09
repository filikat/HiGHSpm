#include "Ipm.h"

#include <cassert>
#include <cmath>
#include <iostream>

// =======================================================================
// LOAD THE PROBLEM
// =======================================================================
void Ipm::load(const int num_var, const int num_con, const double* obj,
               const double* rhs, const double* lower, const double* upper,
               const int* A_colptr, const int* A_rowind, const double* A_values,
               const int* constraints, const std::string& pb_name) {
  if (!obj || !rhs || !lower || !upper || !A_colptr || !A_rowind || !A_values ||
      !constraints)
    return;

  // count how many slacks are needed
  int num_slacks{};
  for (int i = 0; i < num_con; ++i) {
    if (constraints[i] != kConstraintTypeEqual) {
      ++num_slacks;

      if (constraints[i] != kConstraintTypeLower &&
          constraints[i] != kConstraintTypeUpper) {
        std::cerr << "Wrong constraint type\n";
        return;
      }
    }
  }

  // create model with correct size
  model_.resize(num_var + num_slacks, num_con);

  // temporary storage of matrix A
  std::vector<int> temp_colptr(num_var + num_slacks + 1, 0);
  std::vector<int> temp_rowind(A_colptr[num_var] + num_slacks, 0);
  std::vector<double> temp_values(A_colptr[num_var] + num_slacks, 0.0);

  for (int i = 0; i < num_var; ++i) {
    // copy vector c
    model_.obj_[i] = obj[i];

    // copy original lower bound
    model_.lower_[i] = lower[i];

    // copy original upper bound
    model_.upper_[i] = upper[i];

    // copy original column pointers
    temp_colptr[i] = A_colptr[i];
  }
  // copy last column pointer
  temp_colptr[num_var] = A_colptr[num_var];

  int A_nnz{A_colptr[num_var]};

  for (int i = 0; i < A_nnz; ++i) {
    // copy original row indices
    temp_rowind[i] = A_rowind[i];

    // copy original values
    temp_values[i] = A_values[i];
  }

  // index of the current slack
  int slack_ind{num_var};

  for (int i = 0; i < num_con; ++i) {
    // copy vector b
    model_.rhs_[i] = rhs[i];

    // if constraint is inequality, add a slack variable
    if (constraints[i] != kConstraintTypeEqual) {
      // lower/upper bound for new slack
      if (constraints[i] == kConstraintTypeLower) {
        model_.lower_[slack_ind] = -kInf;
        model_.upper_[slack_ind] = 0.0;
      } else {
        model_.lower_[slack_ind] = 0.0;
        model_.upper_[slack_ind] = kInf;
      }

      // add column of identity to A
      temp_colptr[slack_ind + 1] = temp_colptr[slack_ind] + 1;
      temp_rowind[A_nnz] = i;
      temp_values[A_nnz] = 1.0;
      ++A_nnz;

      model_.obj_[slack_ind] = 0.0;

      ++slack_ind;
    }
  }

  model_.A_.num_col_ = num_var + num_slacks;
  model_.A_.num_row_ = num_con;
  model_.A_.start_ = temp_colptr;
  model_.A_.index_ = temp_rowind;
  model_.A_.value_ = temp_values;

  m_ = model_.num_con_;
  n_ = model_.num_var_;

  model_.pb_name_ = pb_name;
  model_ready_ = true;
}

// =======================================================================
// SOLVE THE LP
// =======================================================================
Output Ipm::solve() {
  // solve only if model is loaded
  if (!model_ready_) return Output{};

  printf("-----------------------------------------\n");
  printf("Problem %s\n", model_.pb_name_.c_str());
  printf("\trows %d\n", model_.A_.num_row_);
  printf("\tcols %d\n", model_.A_.num_col_);
  printf("\tnnz  %d\n", model_.A_.numNz());
  printf("\tusing %s\n", option_nla_ == kOptionNlaAugmented
                             ? "augmented system"
                             : "normal equations");
  printf("-----------------------------------------\n");

  double timer_iterations = getWallTime();

  // ------------------------------------------
  // ---- INITIALIZE --------------------------
  // ------------------------------------------

  // iterations counter
  int iter{};

  // initialize stepsize
  double alpha_primal{};
  double alpha_dual{};

  // initialize starting point
  it_ = Iterate(m_, n_);

  // initialize linear solver
  CholmodSolver cholmod_solver;
  MA57Solver ma57_solver;
  MA86Solver ma86_solver;
  MA87Solver ma87_solver;
  MA97Solver ma97_solver;
  HFactorSolver hfactor_solver;
  FactorHiGHSSolver factorHiGHS_solver;
  // linsol_ = &cholmod_solver;
  linsol_ = &factorHiGHS_solver;
  // linsol_ = &hfactor_solver;

  linsol2_ = &cholmod_solver;
  // linsol2_ = &hfactor_solver;

  // perform any preliminary calculations for the linear solver
  linsol_->setup(model_.A_, option_nla_);
  linsol2_->setup(model_.A_, option_nla_);

  // comment the following line to use default starting point
  computeStartingPoint();

  // initialize residuals
  Residuals res(m_, n_);
  computeResiduals1234(res);

  // initialize infeasibilities and mu
  double primal_infeas = norm2(res.res1) / (1 + norm2(model_.rhs_));
  double dual_infeas = norm2(res.res4) / (1 + norm2(model_.obj_));
  double mu = computeMu();
  double primal_obj{};
  double dual_obj{};
  double pd_gap = 1.0;

  std::string status;

  // ------------------------------------------
  // ---- MAIN LOOP ---------------------------
  // ------------------------------------------

  printf("\n");
  while (iter < kMaxIterations) {
    // Check that iterate is not NaN or Inf
    assert(!it_.isNaN());
    assert(!it_.isInf());

    // Stopping criterion
    if (iter > 0 && pd_gap < kIpmTolerance &&  // primal-dual gap is small
        mu < kIpmTolerance &&                  // mu
        primal_infeas < kIpmTolerance &&       // primal feasibility
        dual_infeas < kIpmTolerance) {         // dual feasibility
      printf("\n===== Optimal solution found =====\n\n");
      status = "optimal";
      break;
    }

    // Possibly print header
    if (iter % 20 == 0)
      printf(
          " iter         primal obj            dual obj        pinf       dinf "
          "       mu        alpha_p    alpha_d      p/d rel gap      time\n");

    ++iter;

    // Compute diagonal scaling
    std::vector<double> scaling(n_, 0.0);
    computeScaling(scaling);

    // Clear any existing data in the linear solver now that scaling has changed
    linsol_->clear();
    linsol2_->clear();

    // Initialize Newton direction
    NewtonDir delta(m_, n_);

    if (option_predcor_ == 0) {
      bool is_corrector = false;

      // Heuristic to choose sigma
      double sigma{};
      if (iter == 1) {
        sigma = kSigmaInitial;
      } else {
        sigma = pow(std::max(1.0 - alpha_primal, 1.0 - alpha_dual), 5.0);
        sigma = std::max(sigma, kSigmaMin);
        sigma = std::min(sigma, kSigmaMax);
      }

      // Compute last two residuals with correct value of sigma
      computeResiduals56(sigma * mu, delta, is_corrector, res);

      // Solve Newton system
      if (solveNewtonSystem(model_.A_, scaling, res, is_corrector, delta)) {
        std::cerr << "Error while solving Newton system\n";
        status = "Error";
        break;
      }

      // Compute full Newton direction
      recoverDirection(res, is_corrector, delta);

      // CheckResiduals(Delta, Res);

    } else {
      // *********************************************************************
      // PREDICTOR
      // *********************************************************************
      bool is_corrector = false;

      // Compute last two residuals for predictor
      computeResiduals56(0, delta, is_corrector, res);

      // Solve Newton system for predictor
      if (solveNewtonSystem(model_.A_, scaling, res, is_corrector, delta)) {
        std::cerr << "Error while solving Newton system\n";
        status = "Error";
        break;
      }

      // Compute full Newton direction for predictor
      recoverDirection(res, is_corrector, delta);
      // *********************************************************************

      // *********************************************************************
      // CORRECTOR
      // *********************************************************************
      is_corrector = true;

      // Compute sigma based on the predictor
      double sigma = computeSigmaCorrector(delta, mu);

      // Compute last two residuals for corrector
      computeResiduals56(sigma * mu, delta, is_corrector, res);

      // Initialize corrector direction
      NewtonDir delta_cor(m_, n_);

      // Solve Newton system for corrector
      res.res1.assign(m_, 0.0);
      if (solveNewtonSystem(model_.A_, scaling, res, is_corrector, delta_cor)) {
        std::cerr << "Error while solving Newton system\n";
        status = "Error";
        break;
      }

      // Compute full Newton direction for corrector
      recoverDirection(res, is_corrector, delta_cor);

      // Add corrector to predictor
      vectorAdd(delta.x, delta_cor.x, 1.0);
      vectorAdd(delta.y, delta_cor.y, 1.0);
      vectorAdd(delta.xl, delta_cor.xl, 1.0);
      vectorAdd(delta.xu, delta_cor.xu, 1.0);
      vectorAdd(delta.zl, delta_cor.zl, 1.0);
      vectorAdd(delta.zu, delta_cor.zu, 1.0);
      // *********************************************************************
    }

    // Find step-sizes
    computeStepSizes(delta, alpha_primal, alpha_dual);

    // Make the step
    vectorAdd(it_.x, delta.x, alpha_primal);
    vectorAdd(it_.xl, delta.xl, alpha_primal);
    vectorAdd(it_.xu, delta.xu, alpha_primal);
    vectorAdd(it_.y, delta.y, alpha_dual);
    vectorAdd(it_.zl, delta.zl, alpha_dual);
    vectorAdd(it_.zu, delta.zu, alpha_dual);

    // Compute first four residuals of the new iterate
    computeResiduals1234(res);

    // Compute mu of the new iterate
    mu = computeMu();

    // Print output to screen
    primal_infeas = norm2(res.res1) / (1 + norm2(model_.rhs_));
    dual_infeas = norm2(res.res4) / (1 + norm2(model_.obj_));
    primal_obj = dotProd(it_.x, model_.obj_);

    // compute dual objective
    dual_obj = dotProd(it_.y, model_.rhs_);
    for (int i = 0; i < n_; ++i) {
      if (model_.hasLb(i)) {
        dual_obj += model_.lower_[i] * it_.zl[i];
      }
      if (model_.hasUb(i)) {
        dual_obj -= model_.upper_[i] * it_.zu[i];
      }
    }

    // compute scaled primal-dual gap
    pd_gap = std::fabs(primal_obj - dual_obj) /
             (1 + 0.5 * std::fabs(primal_obj + dual_obj));

    printf(
        "%5d %20.10e %20.10e %10.2e %10.2e %10.2e %10.2f %10.2f %15.2e "
        "%10.1fs\n",
        iter, primal_obj, dual_obj, primal_infeas, dual_infeas, mu,
        alpha_primal, alpha_dual, pd_gap, getWallTime() - timer_iterations);

    // it.print(iter);
    // delta.print(iter);
    // res.print(iter);
  }

  if (status.empty()) status = "max iter";

  // output struct
  Output out{};
  out.it = std::move(it_);
  out.iterations = iter;
  out.primal_infeas = primal_infeas;
  out.dual_infeas = dual_infeas;
  out.mu = mu;
  out.status = status;

  return out;
}

// =======================================================================
// COMPUTE MU
// =======================================================================
double Ipm::computeMu() {
  double mu = 0.0;
  int number_finite_bounds{};
  for (int i = 0; i < n_; ++i) {
    if (model_.hasLb(i)) {
      mu += it_.xl[i] * it_.zl[i];
      ++number_finite_bounds;
    }
    if (model_.hasUb(i)) {
      mu += it_.xu[i] * it_.zu[i];
      ++number_finite_bounds;
    }
  }
  return (mu / number_finite_bounds);
}

// =======================================================================
// COMPUTE RESIDUALS
// =======================================================================
void Ipm::computeResiduals1234(Residuals& res) {
  // res1
  res.res1 = model_.rhs_;
  model_.A_.alphaProductPlusY(-1.0, it_.x, res.res1);

  // res2
  for (int i = 0; i < n_; ++i) {
    if (model_.hasLb(i)) {
      res.res2[i] = model_.lower_[i] - it_.x[i] + it_.xl[i];
    } else {
      res.res2[i] = 0.0;
    }
  }

  // res3
  for (int i = 0; i < n_; ++i) {
    if (model_.hasUb(i)) {
      res.res3[i] = model_.upper_[i] - it_.x[i] - it_.xu[i];
    } else {
      res.res3[i] = 0.0;
    }
  }

  // res4
  res.res4 = model_.obj_;
  model_.A_.alphaProductPlusY(-1.0, it_.y, res.res4, true);
  for (int i = 0; i < n_; ++i) {
    if (model_.hasLb(i)) {
      res.res4[i] -= it_.zl[i];
    }
    if (model_.hasUb(i)) {
      res.res4[i] += it_.zu[i];
    }
  }

  // Check for NaN or Inf
  assert(!res.isNaN());
  assert(!res.isInf());
}

void Ipm::computeResiduals56(const double sigma_mu, const NewtonDir& delta_aff,
                             bool is_corrector, Residuals& res) {
  if (!is_corrector) {
    for (int i = 0; i < n_; ++i) {
      // res5
      if (model_.hasLb(i)) {
        res.res5[i] = sigma_mu - it_.xl[i] * it_.zl[i];
      } else {
        res.res5[i] = 0.0;
      }

      // res6
      if (model_.hasUb(i)) {
        res.res6[i] = sigma_mu - it_.xu[i] * it_.zu[i];
      } else {
        res.res6[i] = 0.0;
      }
    }

  } else {
    for (int i = 0; i < n_; ++i) {
      // res5
      if (model_.hasLb(i)) {
        res.res5[i] = sigma_mu - delta_aff.xl[i] * delta_aff.zl[i];
      } else {
        res.res5[i] = 0.0;
      }

      // res6
      if (model_.hasUb(i)) {
        res.res6[i] = sigma_mu - delta_aff.xu[i] * delta_aff.zu[i];
      } else {
        res.res6[i] = 0.0;
      }
    }
  }
}

std::vector<double> Ipm::computeResiduals7(const Residuals& res,
                                           bool is_corrector) {
  std::vector<double> res7;

  if (!is_corrector) {
    res7 = res.res4;
    for (int i = 0; i < n_; ++i) {
      if (model_.hasLb(i)) {
        res7[i] -= ((res.res5[i] + it_.zl[i] * res.res2[i]) / it_.xl[i]);
      }
      if (model_.hasUb(i)) {
        res7[i] += ((res.res6[i] - it_.zu[i] * res.res3[i]) / it_.xu[i]);
      }
    }

  } else {
    res7.resize(n_, 0.0);
    for (int i = 0; i < n_; ++i) {
      if (model_.hasLb(i)) {
        res7[i] -= (res.res5[i] / it_.xl[i]);
      }
      if (model_.hasUb(i)) {
        res7[i] += (res.res6[i] / it_.xu[i]);
      }
    }
  }

  return res7;
}

std::vector<double> Ipm::computeResiduals8(const HighsSparseMatrix& A,
                                           const std::vector<double>& scaling,
                                           const Residuals& res,
                                           const std::vector<double>& res7,
                                           bool is_corrector) {
  std::vector<double> res8;

  if (is_corrector) {
    res8.resize(m_, 0.0);
  } else {
    res8 = res.res1;
  }

  std::vector<double> temp(res7);

  // temp = Theta * res7
  vectorDivide(temp, scaling);

  // res8 += A * temp
  A.alphaProductPlusY(1.0, temp, res8);

  return res8;
}

// =======================================================================
// COMPUTE SCALING
// =======================================================================
void Ipm::computeScaling(std::vector<double>& scaling) {
  for (int i = 0; i < n_; ++i) {
    if (model_.hasLb(i)) {
      scaling[i] += it_.zl[i] / it_.xl[i];
    }
    if (model_.hasUb(i)) {
      scaling[i] += it_.zu[i] / it_.xu[i];
    }

    // add primal regularization
    scaling[i] += kPrimalRegularization;
  }
}

// =======================================================================
// SOLVE NEWTON SYSTEM
// =======================================================================
int Ipm::solveNewtonSystem(const HighsSparseMatrix& A,
                           const std::vector<double>& scaling,
                           const Residuals& res, bool is_corrector,
                           NewtonDir& delta) {
  // Compute residual 7
  std::vector<double> res7{computeResiduals7(res, is_corrector)};

  // Identify whether the augmented system should be solved directly
  const bool use_direct_augmented = option_nla_ == kOptionNlaAugmented;
  // Identify whether the Newton system should be solved directly
  const bool use_direct_newton = option_nla_ == kOptionNlaNormEq;

  // Augmented system is
  //
  // [-scaling A^T][dx] = [res7]
  // [A        0  ][dy]   [res1]
  //
  std::vector<double> theta;
  scaling2theta(scaling, theta);

  const bool first_call_with_theta = !linsol_->valid_;

  if (use_direct_newton) {
    // Have to solve the Newton system
    // A.scaling.A^T delta.y = res8

    // Compute res8
    std::vector<double> res8{
        computeResiduals8(A, scaling, res, res7, is_corrector)};

    if (first_call_with_theta) {
      int newton_invert_status = linsol_->factorNE(A, theta);
      linsol2_->factorNE(A, theta);

      if (newton_invert_status) return newton_invert_status;
    }

    int newton_solve_status = linsol_->solveNE(A, theta, res8, delta.y);

    if (newton_solve_status) return newton_solve_status;

    {
      std::vector<double> temp(m_, 0.0);
      linsol2_->solveNE(A, theta, res8, temp);

      double diff{};
      double normy{};
      for (int i = 0; i < m_; ++i) {
        diff += (temp[i] - delta.y[i]) * (temp[i] - delta.y[i]);
        normy += delta.y[i] * delta.y[i];
      }
      diff = sqrt(diff);
      normy = sqrt(normy);

      printf("Error: %e\n", diff / normy);
    }

    // Compute delta.x
    // *********************************************************************
    // Deltax = A^T * Deltay - res7;
    delta.x = res7;
    A.alphaProductPlusY(-1.0, delta.y, delta.x, true);
    vectorScale(delta.x, -1.0);

    // Deltax = Theta * Deltax
    vectorDivide(delta.x, scaling);
    // *********************************************************************
  }

  if (use_direct_augmented) {
    if (first_call_with_theta) {
      int augmented_invert_status = linsol_->factorAS(A, theta);
      linsol2_->factorAS(A, theta);
      if (augmented_invert_status) return augmented_invert_status;
    }

    // When solving the augmented system, the right side for the
    // predictor should be (res7; res1), but for the corrector it
    // should be (res7; 0)
    if (is_corrector) {
      // Check that res1 has been zeroed
      const double res1_norm = norm2(res.res1);
      assert(!res1_norm);
    }

    linsol_->solveAS(A, theta, res7, res.res1, delta.x, delta.y);

    {
      std::vector<double> tempx(n_, 0.0);
      std::vector<double> tempy(m_, 0.0);
      linsol2_->solveAS(A, theta, res7, res.res1, tempx, tempy);

      double diff{};
      double normy{};
      for (int i = 0; i < m_; ++i) {
        diff += (tempy[i] - delta.y[i]) * (tempy[i] - delta.y[i]);
        normy += delta.y[i] * delta.y[i];
      }
      for (int i = 0; i < n_; ++i) {
        diff += (tempx[i] - delta.x[i]) * (tempx[i] - delta.x[i]);
        normy += delta.x[i] * delta.x[i];
      }
      diff = sqrt(diff);
      normy = sqrt(normy);

      printf("Error: %e\n", diff / normy);
    }
  }

  return 0;
}

// =======================================================================
// FULL NEWTON DIRECTION
// =======================================================================
void Ipm::recoverDirection(const Residuals& res, bool is_corrector,
                           NewtonDir& delta) {
  if (!is_corrector) {
    // Deltaxl
    delta.xl = delta.x;
    vectorAdd(delta.xl, res.res2, -1.0);

    // Deltaxu
    delta.xu = res.res3;
    vectorAdd(delta.xu, delta.x, -1.0);
  } else {
    // Deltaxl
    delta.xl = delta.x;

    // Deltaxu
    delta.xu = delta.x;
    vectorScale(delta.xu, -1.0);
  }

  // Deltazl
  delta.zl = res.res5;
  vectorAddMult(delta.zl, it_.zl, delta.xl, -1.0);
  vectorDivide(delta.zl, it_.xl);

  // Deltazu
  delta.zu = res.res6;
  vectorAddMult(delta.zu, it_.zu, delta.xu, -1.0);
  vectorDivide(delta.zu, it_.xu);

  // Check for NaN of Inf
  assert(!delta.isNaN());
  assert(!delta.isInf());
}

// =======================================================================
// COMPUTE STEP-SIZES
// =======================================================================
void Ipm::computeStepSizes(const NewtonDir& delta, double& alpha_primal,
                           double& alpha_dual) {
  alpha_primal = 1.0;
  for (int i = 0; i < n_; ++i) {
    if (delta.xl[i] < 0 && model_.hasLb(i)) {
      alpha_primal = std::min(alpha_primal, -it_.xl[i] / delta.xl[i]);
    }
    if (delta.xu[i] < 0 && model_.hasUb(i)) {
      alpha_primal = std::min(alpha_primal, -it_.xu[i] / delta.xu[i]);
    }
  }
  alpha_primal *= kInteriorScaling;

  alpha_dual = 1.0;
  for (int i = 0; i < n_; ++i) {
    if (delta.zl[i] < 0 && model_.hasLb(i)) {
      alpha_dual = std::min(alpha_dual, -it_.zl[i] / delta.zl[i]);
    }
    if (delta.zu[i] < 0 && model_.hasUb(i)) {
      alpha_dual = std::min(alpha_dual, -it_.zu[i] / delta.zu[i]);
    }
  }
  alpha_dual *= kInteriorScaling;

  assert(alpha_primal > 0 && alpha_primal < 1 && alpha_dual > 0 &&
         alpha_dual < 1);
}

// ===================================================================================
// COMPUTE STARTING POINT
// ===================================================================================
void Ipm::computeStartingPoint() {
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
  std::vector<double> temp_m(m_);

  if (option_nla_ == kOptionNlaNormEq) {
    // use y to store b-A*x
    it_.y = model_.rhs_;
    model_.A_.alphaProductPlusY(-1.0, it_.x, it_.y);

    // solve A*A^T * dx = b-A*x with factorization and store the result in
    // temp_m

    // factorize A*A^T
    int newton_invert_status = linsol_->factorNE(model_.A_, temp_scaling);
    if (newton_invert_status) std::cerr << "Error in factorization of A*A^T\n";

    int newton_solve_status =
        linsol_->solveNE(model_.A_, temp_scaling, it_.y, temp_m);
    if (newton_solve_status) std::cerr << "Error in solution of A*A^T\n";

  } else if (option_nla_ == kOptionNlaAugmented) {
    // obtain solution of A*A^T * dx = b-A*x by solving
    // [ -I  A^T] [...] = [ -x]
    // [  A   0 ] [ dx] = [ b ]

    int augmented_invert_status = linsol_->factorAS(model_.A_, temp_scaling);
    if (augmented_invert_status)
      std::cerr << "Error in factorization of A*A^T\n";

    std::vector<double> rhs_x(n_);
    for (int i = 0; i < n_; ++i) rhs_x[i] = -it_.x[i];
    std::vector<double> lhs_x(n_);
    int augmented_solve_status = linsol_->solveAS(
        model_.A_, temp_scaling, rhs_x, model_.rhs_, lhs_x, temp_m);

    if (augmented_solve_status) std::cerr << "Error in solution of A*A^T\n";
  }

  // compute dx = A^T * (A*A^T)^{-1} * (b-A*x) and store the result in xl
  std::fill(it_.xl.begin(), it_.xl.end(), 0.0);
  model_.A_.alphaProductPlusY(1.0, temp_m, it_.xl, true);

  // x += dx;
  vectorAdd(it_.x, it_.xl, 1.0);
  // *********************************************************************

  // *********************************************************************
  // xl, xu starting point
  // *********************************************************************
  // compute xl, xu that satisfy linear constraints
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
  // *********************************************************************

  // *********************************************************************
  // y starting point
  // *********************************************************************

  if (option_nla_ == kOptionNlaNormEq) {
    // compute A*c
    std::fill(temp_m.begin(), temp_m.end(), 0.0);
    model_.A_.alphaProductPlusY(1.0, model_.obj_, temp_m);

    int newton_solve_status =
        linsol_->solveNE(model_.A_, temp_scaling, temp_m, it_.y);
    if (newton_solve_status) std::cerr << "Error in solution of A*A^T\n";

  } else if (option_nla_ == kOptionNlaAugmented) {
    // obtain solution of A*A^T * y = A*c by solving
    // [ -I  A^T] [...] = [ c ]
    // [  A   0 ] [ y ] = [ 0 ]

    std::vector<double> rhs_y(m_, 0.0);
    std::vector<double> lhs_x(n_);

    int augmented_solve_status = linsol_->solveAS(
        model_.A_, temp_scaling, model_.obj_, rhs_y, lhs_x, it_.y);
    if (augmented_solve_status) std::cerr << "Error in solution of A*A^T\n";
  }
  // *********************************************************************

  // *********************************************************************
  // zl, zu starting point
  // *********************************************************************
  // compute c - A^T * y and store in zl
  it_.zl = model_.obj_;
  model_.A_.alphaProductPlusY(-1.0, it_.y, it_.zl, true);

  // split result between zl and zu
  violation = 0.0;
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
  // *********************************************************************

  // *********************************************************************
  // improve centrality
  // *********************************************************************
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
  // *********************************************************************
}

double Ipm::computeSigmaCorrector(const NewtonDir& delta_aff, double mu) {
  // stepsizes of predictor direction
  double alpha_p{};
  double alpha_d{};
  computeStepSizes(delta_aff, alpha_p, alpha_d);

  // mu using predictor direction
  double mu_aff = 0.0;
  int number_finite_bounds{};
  for (int i = 0; i < n_; ++i) {
    if (model_.hasLb(i)) {
      mu_aff += (it_.xl[i] + alpha_p * delta_aff.xl[i]) *
                (it_.zl[i] + alpha_d * delta_aff.zl[i]);
      ++number_finite_bounds;
    }
    if (model_.hasUb(i)) {
      mu_aff += (it_.xu[i] + alpha_p * delta_aff.xu[i]) *
                (it_.zu[i] + alpha_d * delta_aff.zu[i]);
      ++number_finite_bounds;
    }
  }
  mu_aff /= number_finite_bounds;

  // heuristic to choose sigma
  double ratio = mu_aff / mu;

  return ratio * ratio * ratio;
}

void Ipm::checkResiduals(const NewtonDir& delta, const Residuals& res) const {
  std::vector<double> temp_m(m_, 0.0);
  std::vector<double> temp_n(n_, 0.0);

  // check A * delta.x - res1
  model_.A_.product(temp_m, delta.x);
  double norm_diff_1 = infNormDiff(temp_m, res.res1);
  vectorAdd(temp_m, res.res1, -1.0);
  vectorDivide(temp_m, res.res1);
  double rel_norm_1 = infNorm(temp_m);

  // check delta.x - delta.xl - res2
  temp_n = delta.x;
  vectorAdd(temp_n, delta.xl, -1.0);
  double norm_diff_2 = infNormDiff(temp_n, res.res2);
  vectorAdd(temp_n, res.res2, -1.0);
  vectorDivide(temp_n, res.res2);
  double rel_norm_2 = infNorm(temp_n);

  // check delta.x + delta.xu - res3
  temp_n = delta.x;
  vectorAdd(temp_n, delta.xu);
  double norm_diff_3 = infNormDiff(temp_n, res.res3);
  vectorAdd(temp_n, res.res3, -1.0);
  vectorDivide(temp_n, res.res3);
  double rel_norm_3 = infNorm(temp_n);

  // check A^T * delta.y + delta.zl - delta.zu - res4
  std::fill(temp_n.begin(), temp_n.end(), 0.0);
  model_.A_.productTranspose(temp_n, delta.y);
  vectorAdd(temp_n, delta.zl);
  vectorAdd(temp_n, delta.zu, -1.0);
  double norm_diff_4 = infNormDiff(temp_n, res.res4);
  vectorAdd(temp_n, res.res4, -1.0);
  vectorDivide(temp_n, res.res4);
  double rel_norm_4 = infNorm(temp_n);

  // check Zl * delta.xl + Xl * delta.zl - res5
  temp_n = delta.xl;
  vectorMultiply(temp_n, it_.zl);
  vectorAddMult(temp_n, delta.zl, it_.xl);
  double norm_diff_5 = infNormDiff(temp_n, res.res5);
  vectorAdd(temp_n, res.res5, -1.0);
  vectorDivide(temp_n, res.res5);
  double rel_norm_5 = infNorm(temp_n);

  // check Zu * delta.xu + Xu * delta.zu - res6
  temp_n = delta.xu;
  vectorMultiply(temp_n, it_.zu);
  vectorAddMult(temp_n, delta.zu, it_.xu);
  double norm_diff_6 = infNormDiff(temp_n, res.res6);
  vectorAdd(temp_n, res.res6, -1.0);
  vectorDivide(temp_n, res.res6);
  double rel_norm_6 = infNorm(temp_n);

  printf("Diff norms:\n");
  printf("1: %.2e\n", norm_diff_1);
  printf("2: %.2e\n", norm_diff_2);
  printf("3: %.2e\n", norm_diff_3);
  printf("4: %.2e\n", norm_diff_4);
  printf("5: %.2e\n", norm_diff_5);
  printf("6: %.2e\n", norm_diff_6);
}
