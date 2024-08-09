#include "IPM_caller.h"

#include <cassert>
#include <cmath>
#include <iostream>

// =======================================================================
// LOAD THE PROBLEM
// =======================================================================
void IPM_caller::Load(const int num_var, const int num_con, const double *obj,
                      const double *rhs, const double *lower,
                      const double *upper, const int *A_colptr,
                      const int *A_rowind, const double *A_values,
                      const int *constraints, const std::string &pb_name) {
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
  model.resize(num_var + num_slacks, num_con);

  // temporary storage of matrix A
  std::vector<int> temp_colptr(num_var + num_slacks + 1, 0);
  std::vector<int> temp_rowind(A_colptr[num_var] + num_slacks, 0);
  std::vector<double> temp_values(A_colptr[num_var] + num_slacks, 0.0);

  for (int i = 0; i < num_var; ++i) {
    // copy vector c
    model.obj[i] = obj[i];

    // copy original lower bound
    model.lower[i] = lower[i];

    // copy original upper bound
    model.upper[i] = upper[i];

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
    model.rhs[i] = rhs[i];

    // if constraint is inequality, add a slack variable
    if (constraints[i] != kConstraintTypeEqual) {
      // lower/upper bound for new slack
      if (constraints[i] == kConstraintTypeLower) {
        model.lower[slack_ind] = -INF;
        model.upper[slack_ind] = 0.0;
      } else {
        model.lower[slack_ind] = 0.0;
        model.upper[slack_ind] = INF;
      }

      // add column of identity to A
      temp_colptr[slack_ind + 1] = temp_colptr[slack_ind] + 1;
      temp_rowind[A_nnz] = i;
      temp_values[A_nnz] = 1.0;
      ++A_nnz;

      model.obj[slack_ind] = 0.0;

      ++slack_ind;
    }
  }

  model.highs_a.num_col_ = num_var + num_slacks;
  model.highs_a.num_row_ = num_con;
  model.highs_a.start_ = temp_colptr;
  model.highs_a.index_ = temp_rowind;
  model.highs_a.value_ = temp_values;
  //  assert(equalMatrix("In load"));

  m = model.num_con;
  n = model.num_var;

  model.pb_name = pb_name;
  model_ready = true;
}

// =======================================================================
// SOLVE THE LP
// =======================================================================
Output IPM_caller::Solve() {
  // solve only if model is loaded
  if (!model_ready)
    return Output{};

  printf("-----------------------------------------\n");
  printf("Problem %s\n", model.pb_name.c_str());
  printf("A rows %d\n", model.highs_a.num_row_);
  printf("A cols %d\n", model.highs_a.num_col_);
  printf("A nnz  %d\n", model.highs_a.numNz());
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
  It = Iterate(m, n);

  // initialize linear solver
  CholmodSolver cholmodsolver;
  MA86Solver ma86solver;
  MA87Solver ma87solver;
  MA97Solver ma97solver;
  //linsol = &cholmodsolver;
  linsol =&ma86solver;

  // comment the following line to use default starting point
  ComputeStartingPoint();

  // initialize residuals
  Residuals Res(m, n);
  ComputeResiduals_1234(Res);

  // initialize infeasibilities and mu
  double primal_infeas = Norm2(Res.res1) / (1 + Norm2(model.rhs));
  double dual_infeas = Norm2(Res.res4) / (1 + Norm2(model.obj));
  double mu = ComputeMu();
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
    assert(!It.isNaN());
    assert(!It.isInf());

    // Stopping criterion
    if (iter > 0 && pd_gap < kIpmTolerance && // primal-dual gap is small
        mu < kIpmTolerance &&                 // mu
        primal_infeas < kIpmTolerance &&      // primal feasibility
        dual_infeas < kIpmTolerance) {        // dual feasibility
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
    std::vector<double> scaling(n, 0.0);
    ComputeScaling(scaling);

    // Clear any existing INVERT now that scaling has changed
    linsol->Clear();

    // Initialize Newton direction
    NewtonDir Delta(m, n);

    if (option_predcor == 0) {
      bool isCorrector = false;

      // Heuristic to choose sigma
      double sigma{};
      if (iter == 1) {
        sigma = sigma_i;
      } else {
        sigma = pow(std::max(1.0 - alpha_primal, 1.0 - alpha_dual), 5.0);
        sigma = std::max(sigma, sigma_min);
        sigma = std::min(sigma, sigma_max);
      }

      // Compute last two residuals with correct value of sigma
      ComputeResiduals_56(sigma * mu, Delta, isCorrector, Res);

      // Solve Newton system
      if (SolveNewtonSystem(model.highs_a, scaling, Res, isCorrector, Delta)) {
        std::cerr << "Error while solving Newton system\n";
        status = "Error";
        break;
      }

      // Compute full Newton direction
      RecoverDirection(Res, isCorrector, Delta);

      // CheckResiduals(Delta, Res);

    } else {
      // *********************************************************************
      // PREDICTOR
      // *********************************************************************
      bool isCorrector = false;

      // Compute last two residuals for predictor
      ComputeResiduals_56(0, Delta, isCorrector, Res);

      // Solve Newton system for predictor
      if (SolveNewtonSystem(model.highs_a, scaling, Res, isCorrector, Delta)) {
        std::cerr << "Error while solving Newton system\n";
        status = "Error";
        break;
      }

      // Compute full Newton direction for predictor
      RecoverDirection(Res, isCorrector, Delta);
      // *********************************************************************

      // *********************************************************************
      // CORRECTOR
      // *********************************************************************
      isCorrector = true;

      // Compute sigma based on the predictor
      double sigma = ComputeSigmaCorrector(Delta, mu);

      // Compute last two residuals for corrector
      ComputeResiduals_56(sigma * mu, Delta, isCorrector, Res);

      // Initialize corrector direction
      NewtonDir DeltaCor(m, n);

      // Solve Newton system for corrector
      Res.res1.assign(m, 0.0);
      if (SolveNewtonSystem(model.highs_a, scaling, Res, isCorrector,
                            DeltaCor)) {
        std::cerr << "Error while solving Newton system\n";
        status = "Error";
        break;
      }

      // Compute full Newton direction for corrector
      RecoverDirection(Res, isCorrector, DeltaCor);

      // Add corrector to predictor
      // (all these add can be avoided with a different implementation, but for
      // now it works)
      VectorAdd(Delta.x, DeltaCor.x, 1.0);
      VectorAdd(Delta.y, DeltaCor.y, 1.0);
      VectorAdd(Delta.xl, DeltaCor.xl, 1.0);
      VectorAdd(Delta.xu, DeltaCor.xu, 1.0);
      VectorAdd(Delta.zl, DeltaCor.zl, 1.0);
      VectorAdd(Delta.zu, DeltaCor.zu, 1.0);
      // *********************************************************************
    }

    // Find step-sizes
    ComputeStepSizes(Delta, alpha_primal, alpha_dual);

    // Make the step
    VectorAdd(It.x, Delta.x, alpha_primal);
    VectorAdd(It.xl, Delta.xl, alpha_primal);
    VectorAdd(It.xu, Delta.xu, alpha_primal);
    VectorAdd(It.y, Delta.y, alpha_dual);
    VectorAdd(It.zl, Delta.zl, alpha_dual);
    VectorAdd(It.zu, Delta.zu, alpha_dual);

    // Compute first four residuals of the new iterate
    ComputeResiduals_1234(Res);

    // Compute mu of the new iterate
    mu = ComputeMu();

    // Print output to screen
    primal_infeas = Norm2(Res.res1) / (1 + Norm2(model.rhs));
    dual_infeas = Norm2(Res.res4) / (1 + Norm2(model.obj));
    primal_obj = DotProd(It.x, model.obj);

    // compute dual objective
    dual_obj = DotProd(It.y, model.rhs);
    for (int i = 0; i < n; ++i) {
      if (model.has_lb(i)) {
        dual_obj += model.lower[i] * It.zl[i];
      }
      if (model.has_ub(i)) {
        dual_obj -= model.upper[i] * It.zu[i];
      }
    }

    // compute scaled primal-dual gap
    pd_gap = std::fabs(primal_obj - dual_obj) /
             (1 + 0.5 * std::fabs(primal_obj + dual_obj));

    printf("%5d %20.10e %20.10e %10.2e %10.2e %10.2e %10.2f %10.2f %15.2e "
           "%10.1fs\n",
           iter, primal_obj, dual_obj, primal_infeas, dual_infeas, mu,
           alpha_primal, alpha_dual, pd_gap, getWallTime() - timer_iterations);

    // It.print(iter);
    // Delta.print(iter);
    // Res.print(iter);
  }

  if (status.empty())
    status = "max iter";

  // output struct
  Output out{};
  out.It = It;
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
double IPM_caller::ComputeMu() {
  double mu = 0.0;
  int number_finite_bounds{};
  for (int i = 0; i < n; ++i) {
    if (model.has_lb(i)) {
      mu += It.xl[i] * It.zl[i];
      ++number_finite_bounds;
    }
    if (model.has_ub(i)) {
      mu += It.xu[i] * It.zu[i];
      ++number_finite_bounds;
    }
  }
  return (mu / number_finite_bounds);
}

// =======================================================================
// COMPUTE RESIDUALS
// =======================================================================
void IPM_caller::ComputeResiduals_1234(Residuals &Res) {
  // res1
  Res.res1 = model.rhs;
  model.highs_a.alphaProductPlusY(-1.0, It.x, Res.res1);

  // res2
  for (int i = 0; i < n; ++i) {
    if (model.has_lb(i)) {
      Res.res2[i] = model.lower[i] - It.x[i] + It.xl[i];
    } else {
      Res.res2[i] = 0.0;
    }
  }

  // res3
  for (int i = 0; i < n; ++i) {
    if (model.has_ub(i)) {
      Res.res3[i] = model.upper[i] - It.x[i] - It.xu[i];
    } else {
      Res.res3[i] = 0.0;
    }
  }

  // res4
  Res.res4 = model.obj;
  model.highs_a.alphaProductPlusY(-1.0, It.y, Res.res4, true);
  for (int i = 0; i < n; ++i) {
    if (model.has_lb(i)) {
      Res.res4[i] -= It.zl[i];
    }
    if (model.has_ub(i)) {
      Res.res4[i] += It.zu[i];
    }
  }

  // Check for NaN or Inf
  assert(!Res.isNaN());
  assert(!Res.isInf());
}

void IPM_caller::ComputeResiduals_56(const double sigmaMu,
                                     const NewtonDir &DeltaAff,
                                     bool isCorrector, Residuals &Res) {
  if (!isCorrector) {
    for (int i = 0; i < n; ++i) {
      // res5
      if (model.has_lb(i)) {
        Res.res5[i] = sigmaMu - It.xl[i] * It.zl[i];
      } else {
        Res.res5[i] = 0.0;
      }

      // res6
      if (model.has_ub(i)) {
        Res.res6[i] = sigmaMu - It.xu[i] * It.zu[i];
      } else {
        Res.res6[i] = 0.0;
      }
    }

  } else {
    for (int i = 0; i < n; ++i) {
      // res5
      if (model.has_lb(i)) {
        Res.res5[i] = sigmaMu - DeltaAff.xl[i] * DeltaAff.zl[i];
      } else {
        Res.res5[i] = 0.0;
      }

      // res6
      if (model.has_ub(i)) {
        Res.res6[i] = sigmaMu - DeltaAff.xu[i] * DeltaAff.zu[i];
      } else {
        Res.res6[i] = 0.0;
      }
    }
  }
}

std::vector<double> IPM_caller::ComputeResiduals_7(const Residuals &Res,
                                                   bool isCorrector) {
  std::vector<double> res7;

  if (!isCorrector) {
    res7 = Res.res4;
    for (int i = 0; i < n; ++i) {
      if (model.has_lb(i)) {
        res7[i] -= ((Res.res5[i] + It.zl[i] * Res.res2[i]) / It.xl[i]);
      }
      if (model.has_ub(i)) {
        res7[i] += ((Res.res6[i] - It.zu[i] * Res.res3[i]) / It.xu[i]);
      }
    }

  } else {
    res7.resize(n, 0.0);
    for (int i = 0; i < n; ++i) {
      if (model.has_lb(i)) {
        res7[i] -= (Res.res5[i] / It.xl[i]);
      }
      if (model.has_ub(i)) {
        res7[i] += (Res.res6[i] / It.xu[i]);
      }
    }
  }

  return res7;
}

std::vector<double> IPM_caller::ComputeResiduals_8(
    const HighsSparseMatrix &highs_a, const std::vector<double> &scaling,
    const Residuals &Res, const std::vector<double> &res7, bool isCorrector) {
  std::vector<double> res8;

  if (isCorrector) {
    res8.resize(m, 0.0);
  } else {
    res8 = Res.res1;
  }

  std::vector<double> temp(res7);

  // temp = Theta * res7
  VectorDivide(temp, scaling);

  // res8 += A * temp
  highs_a.alphaProductPlusY(1.0, temp, res8);

  return res8;
}

// =======================================================================
// COMPUTE SCALING
// =======================================================================
void IPM_caller::ComputeScaling(std::vector<double> &scaling) {
  for (int i = 0; i < n; ++i) {
    if (model.has_lb(i)) {
      scaling[i] += It.zl[i] / It.xl[i];
    }
    if (model.has_ub(i)) {
      scaling[i] += It.zu[i] / It.xu[i];
    }

    // add primal regularization
    scaling[i] += kPrimalRegularization;
  }
}

// =======================================================================
// SOLVE NEWTON SYSTEM
// =======================================================================
int IPM_caller::SolveNewtonSystem(const HighsSparseMatrix &highs_a,
                                  const std::vector<double> &scaling,
                                  const Residuals &Res, bool isCorrector,
                                  NewtonDir &Delta) {
  // Compute residual 7
  std::vector<double> res7{ComputeResiduals_7(Res, isCorrector)};

  // Identify whether the augmented system should be solved directly
  const bool use_direct_augmented = option_nla == kOptionNlaAugmented;
  // Identify whether the Newton system should be solved directly
  const bool use_direct_newton = option_nla == kOptionNlaNewton;

  // Augmented system is
  //
  // [-scaling A^T][dx] = [res7]
  // [A        0  ][dy]   [res1]
  //
  std::vector<double> theta;
  scaling2theta(scaling, theta);

  const bool first_call_with_theta = !linsol->valid;

  if (use_direct_newton) {
    // Have to solve the Newton system
    // A.scaling.A^T Delta.y = res8

    // Compute res8
    std::vector<double> res8{
        ComputeResiduals_8(highs_a, scaling, Res, res7, isCorrector)};

    if (first_call_with_theta) {
      int newton_invert_status = linsol->FactorNE(highs_a, theta);

      if (newton_invert_status)
        return newton_invert_status;
    }

    int newton_solve_status = linsol->SolveNE(highs_a, theta, res8, Delta.y);

    if (newton_solve_status)
      return newton_solve_status;

    // Compute Delta.x
    // *********************************************************************
    // Deltax = A^T * Deltay - res7;
    Delta.x = res7;
    highs_a.alphaProductPlusY(-1.0, Delta.y, Delta.x, true);
    VectorScale(Delta.x, -1.0);

    // Deltax = Theta * Deltax
    VectorDivide(Delta.x, scaling);
    // *********************************************************************
  }

  if (use_direct_augmented) {

    if (first_call_with_theta) {
      int augmented_invert_status = linsol->FactorAS(highs_a, theta);
      if (augmented_invert_status)
        return augmented_invert_status;
    }

    // When solving the augmented system, the right side for the
    // predictor should be (res7; res1), but for the corrector it
    // should be (res7; 0)
    if (isCorrector) {
      // Check that res1 has been zeroed
      const double res1_norm = Norm2(Res.res1);
      assert(!res1_norm);
    }

    linsol->SolveAS(highs_a, theta, res7, Res.res1, Delta.x, Delta.y);
  }

  return 0;
}

// =======================================================================
// FULL NEWTON DIRECTION
// =======================================================================
void IPM_caller::RecoverDirection(const Residuals &Res, bool isCorrector,
                                  NewtonDir &Delta) {
  if (!isCorrector) {
    // Deltaxl
    Delta.xl = Delta.x;
    VectorAdd(Delta.xl, Res.res2, -1.0);

    // Deltaxu
    Delta.xu = Res.res3;
    VectorAdd(Delta.xu, Delta.x, -1.0);
  } else {
    // Deltaxl
    Delta.xl = Delta.x;

    // Deltaxu
    Delta.xu = Delta.x;
    VectorScale(Delta.xu, -1.0);
  }

  // Deltazl
  Delta.zl = Res.res5;
  VectorAddMult(Delta.zl, It.zl, Delta.xl, -1.0);
  VectorDivide(Delta.zl, It.xl);

  // Deltazu
  Delta.zu = Res.res6;
  VectorAddMult(Delta.zu, It.zu, Delta.xu, -1.0);
  VectorDivide(Delta.zu, It.xu);

  // Check for NaN of Inf
  assert(!Delta.isNaN());
  assert(!Delta.isInf());
}

// =======================================================================
// COMPUTE STEP-SIZES
// =======================================================================
void IPM_caller::ComputeStepSizes(const NewtonDir &Delta, double &alpha_primal,
                                  double &alpha_dual) {
  alpha_primal = 1.0;
  for (int i = 0; i < n; ++i) {
    if (Delta.xl[i] < 0 && model.has_lb(i)) {
      alpha_primal = std::min(alpha_primal, -It.xl[i] / Delta.xl[i]);
    }
    if (Delta.xu[i] < 0 && model.has_ub(i)) {
      alpha_primal = std::min(alpha_primal, -It.xu[i] / Delta.xu[i]);
    }
  }
  alpha_primal *= alpha_interior_scaling;

  alpha_dual = 1.0;
  for (int i = 0; i < n; ++i) {
    if (Delta.zl[i] < 0 && model.has_lb(i)) {
      alpha_dual = std::min(alpha_dual, -It.zl[i] / Delta.zl[i]);
    }
    if (Delta.zu[i] < 0 && model.has_ub(i)) {
      alpha_dual = std::min(alpha_dual, -It.zu[i] / Delta.zu[i]);
    }
  }
  alpha_dual *= alpha_interior_scaling;

  assert(alpha_primal > 0 && alpha_primal < 1 && alpha_dual > 0 &&
         alpha_dual < 1);
}

// ===================================================================================
// COMPUTE STARTING POINT
// ===================================================================================
void IPM_caller::ComputeStartingPoint() {

  // *********************************************************************
  // x starting point
  // *********************************************************************
  // compute feasible x
  for (int i = 0; i < n; ++i) {
    It.x[i] = 0.0;
    It.x[i] = std::max(It.x[i], model.lower[i]);
    It.x[i] = std::min(It.x[i], model.upper[i]);
  }

  const std::vector<double> temp_scaling(n, 1.0);
  std::vector<double> temp_m(m);

  if (option_nla == kOptionNlaNewton) {
    // use y to store b-A*x
    It.y = model.rhs;
    model.highs_a.alphaProductPlusY(-1.0, It.x, It.y);

    // solve A*A^T * dx = b-A*x with factorization and store the result in
    // temp_m

    // factorize A*A^T

    int newton_invert_status = linsol->FactorNE(model.highs_a, temp_scaling);
    if (newton_invert_status)
      std::cerr << "Error in factorization of A*A^T\n";

    int newton_solve_status =
        linsol->SolveNE(model.highs_a, temp_scaling, It.y, temp_m);
    if (newton_solve_status)
      std::cerr << "Error in solution of A*A^T\n";

  } else if (option_nla == kOptionNlaAugmented) {
    // obtain solution of A*A^T * dx = b-A*x by solving
    // [ -I  A^T] [...] = [ -x]
    // [  A   0 ] [ dx] = [ b ]

    int augmented_invert_status = linsol->FactorAS(model.highs_a, temp_scaling);
    if (augmented_invert_status)
      std::cerr << "Error in factorization of A*A^T\n";

    std::vector<double> rhs_x(n);
    for (int i = 0; i < n; ++i)
      rhs_x[i] = -It.x[i];
    std::vector<double> lhs_x(n);
    int augmented_solve_status = linsol->SolveAS(
        model.highs_a, temp_scaling, rhs_x, model.rhs, lhs_x, temp_m);

    if (augmented_solve_status)
      std::cerr << "Error in solution of A*A^T\n";
  }

  // compute dx = A^T * (A*A^T)^{-1} * (b-A*x) and store the result in xl
  std::fill(It.xl.begin(), It.xl.end(), 0.0);
  model.highs_a.alphaProductPlusY(1.0, temp_m, It.xl, true);

  // x += dx;
  VectorAdd(It.x, It.xl, 1.0);
  // *********************************************************************

  // *********************************************************************
  // xl, xu starting point
  // *********************************************************************
  // compute xl, xu that satisfy linear constraints
  double violation{};
  for (int i = 0; i < n; ++i) {
    if (model.has_lb(i)) {
      It.xl[i] = It.x[i] - model.lower[i];
      violation = std::min(violation, It.xl[i]);
    } else {
      It.xl[i] = 0.0;
    }
    if (model.has_ub(i)) {
      It.xu[i] = model.upper[i] - It.x[i];
      violation = std::min(violation, It.xu[i]);
    } else {
      It.xu[i] = 0.0;
    }
  }

  // shift to be positive
  violation = 1.0 + std::max(0.0, -1.5 * violation);
  VectorAdd(It.xl, violation);
  VectorAdd(It.xu, violation);
  // *********************************************************************

  // *********************************************************************
  // y starting point
  // *********************************************************************

  if (option_nla == kOptionNlaNewton) {
    // compute A*c
    std::fill(temp_m.begin(), temp_m.end(), 0.0);
    model.highs_a.alphaProductPlusY(1.0, model.obj, temp_m);

    int newton_solve_status =
        linsol->SolveNE(model.highs_a, temp_scaling, temp_m, It.y);
    if (newton_solve_status)
      std::cerr << "Error in solution of A*A^T\n";

  } else if (option_nla == kOptionNlaAugmented) {
    // obtain solution of A*A^T * y = A*c by solving
    // [ -I  A^T] [...] = [ c ]
    // [  A   0 ] [ y ] = [ 0 ]

    std::vector<double> rhs_y(m, 0.0);
    std::vector<double> lhs_x(n);

    int augmented_solve_status = linsol->SolveAS(model.highs_a, temp_scaling,
                                                 model.obj, rhs_y, lhs_x, It.y);
    if (augmented_solve_status)
      std::cerr << "Error in solution of A*A^T\n";
  }
  // *********************************************************************

  // *********************************************************************
  // zl, zu starting point
  // *********************************************************************
  // compute c - A^T * y and store in zl
  It.zl = model.obj;
  model.highs_a.alphaProductPlusY(-1.0, It.y, It.zl, true);

  // split result between zl and zu
  violation = 0.0;
  for (int i = 0; i < n; ++i) {
    double val = It.zl[i];
    It.zl[i] = 0.0;
    It.zu[i] = 0.0;

    if (model.has_lb(i) && model.has_ub(i)) {
      It.zl[i] = 0.5 * val;
      It.zu[i] = -0.5 * val;
    } else if (model.has_lb(i)) {
      It.zl[i] = val;
    } else if (model.has_ub(i)) {
      It.zu[i] = -val;
    }

    violation = std::min(violation, It.zl[i]);
    violation = std::min(violation, It.zu[i]);
  }

  // shift to be positive
  violation = 1.0 + std::max(0.0, -1.5 * violation);
  for (int i = 0; i < n; ++i) {
    if (model.has_lb(i)) {
      It.zl[i] += violation;
    }
    if (model.has_ub(i)) {
      It.zu[i] += violation;
    }
  }
  // *********************************************************************

  // *********************************************************************
  // improve centrality
  // *********************************************************************
  double xsum{1.0};
  double zsum{1.0};
  double mu{1.0};

  for (int i = 0; i < n; ++i) {
    if (model.has_lb(i)) {
      xsum += It.xl[i];
      zsum += It.zl[i];
      mu += It.xl[i] * It.zl[i];
    }
    if (model.has_ub(i)) {
      xsum += It.xu[i];
      zsum += It.zu[i];
      mu += It.xu[i] * It.zu[i];
    }
  }

  double dx = 0.5 * mu / zsum;
  double dz = 0.5 * mu / xsum;

  VectorAdd(It.xl, dx);
  VectorAdd(It.xu, dx);
  for (int i = 0; i < n; ++i) {
    if (model.has_lb(i)) {
      It.zl[i] += dz;
    }
    if (model.has_ub(i)) {
      It.zu[i] += dz;
    }
  }
  // *********************************************************************
}

double IPM_caller::ComputeSigmaCorrector(const NewtonDir &DeltaAff, double mu) {
  // stepsizes of predictor direction
  double alpha_p{};
  double alpha_d{};
  ComputeStepSizes(DeltaAff, alpha_p, alpha_d);

  // mu using predictor direction
  double mu_aff = 0.0;
  int number_finite_bounds{};
  for (int i = 0; i < n; ++i) {
    if (model.has_lb(i)) {
      mu_aff += (It.xl[i] + alpha_p * DeltaAff.xl[i]) *
                (It.zl[i] + alpha_d * DeltaAff.zl[i]);
      ++number_finite_bounds;
    }
    if (model.has_ub(i)) {
      mu_aff += (It.xu[i] + alpha_p * DeltaAff.xu[i]) *
                (It.zu[i] + alpha_d * DeltaAff.zu[i]);
      ++number_finite_bounds;
    }
  }
  mu_aff /= number_finite_bounds;

  // heuristic to choose sigma
  double ratio = mu_aff / mu;

  return ratio * ratio * ratio;
}

void IPM_caller::CheckResiduals(const NewtonDir &Delta,
                                const Residuals &Res) const {
  std::vector<double> temp_m(m, 0.0);
  std::vector<double> temp_n(n, 0.0);

  // check A * Delta.x - res1
  model.highs_a.product(temp_m, Delta.x);
  double norm_diff_1 = infNormDiff(temp_m, Res.res1);
  VectorAdd(temp_m, Res.res1, -1.0);
  VectorDivide(temp_m, Res.res1);
  double rel_norm_1 = infNorm(temp_m);

  // check Delta.x - Delta.xl - res2
  temp_n = Delta.x;
  VectorAdd(temp_n, Delta.xl, -1.0);
  double norm_diff_2 = infNormDiff(temp_n, Res.res2);
  VectorAdd(temp_n, Res.res2, -1.0);
  VectorDivide(temp_n, Res.res2);
  double rel_norm_2 = infNorm(temp_n);

  // check Delta.x + Delta.xu - res3
  temp_n = Delta.x;
  VectorAdd(temp_n, Delta.xu);
  double norm_diff_3 = infNormDiff(temp_n, Res.res3);
  VectorAdd(temp_n, Res.res3, -1.0);
  VectorDivide(temp_n, Res.res3);
  double rel_norm_3 = infNorm(temp_n);

  // check A^T * Delta.y + Delta.zl - Delta.zu - res4
  std::fill(temp_n.begin(), temp_n.end(), 0.0);
  model.highs_a.productTranspose(temp_n, Delta.y);
  VectorAdd(temp_n, Delta.zl);
  VectorAdd(temp_n, Delta.zu, -1.0);
  double norm_diff_4 = infNormDiff(temp_n, Res.res4);
  VectorAdd(temp_n, Res.res4, -1.0);
  VectorDivide(temp_n, Res.res4);
  double rel_norm_4 = infNorm(temp_n);

  // check Zl * Delta.xl + Xl * Delta.zl - res5
  temp_n = Delta.xl;
  VectorMultiply(temp_n, It.zl);
  VectorAddMult(temp_n, Delta.zl, It.xl);
  double norm_diff_5 = infNormDiff(temp_n, Res.res5);
  VectorAdd(temp_n, Res.res5, -1.0);
  VectorDivide(temp_n, Res.res5);
  double rel_norm_5 = infNorm(temp_n);

  // check Zu * Delta.xu + Xu * Delta.zu - res6
  temp_n = Delta.xu;
  VectorMultiply(temp_n, It.zu);
  VectorAddMult(temp_n, Delta.zu, It.xu);
  double norm_diff_6 = infNormDiff(temp_n, Res.res6);
  VectorAdd(temp_n, Res.res6, -1.0);
  VectorDivide(temp_n, Res.res6);
  double rel_norm_6 = infNorm(temp_n);

  printf("Diff norms:\n");
  printf("1: %.2e\n", norm_diff_1);
  printf("2: %.2e\n", norm_diff_2);
  printf("3: %.2e\n", norm_diff_3);
  printf("4: %.2e\n", norm_diff_4);
  printf("5: %.2e\n", norm_diff_5);
  printf("6: %.2e\n", norm_diff_6);
}
