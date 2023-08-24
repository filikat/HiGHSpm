#include "IPM_caller.h"
#include <cassert>
#include <cmath>
#include <iostream>

void scaling2theta(const std::vector<double> &scaling,
                   std::vector<double> &theta) {
  const int dim = scaling.size();
  theta.resize(dim);
  for (int i = 0; i < dim; i++)
    theta[i] = 1 / scaling[i];
}

// =======================================================================
// LOAD THE PROBLEM
// =======================================================================
void IPM_caller::Load(const int num_var, const int num_con, const double *obj,
                      const double *rhs, const double *lower,
                      const double *upper, const int *A_colptr,
                      const int *A_rowind, const double *A_values,
                      const int *constraints) {

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
  model_ready = true;
}

// =======================================================================
// SOLVE THE LP
// =======================================================================
Output IPM_caller::Solve() {

  // solve only if model is loaded
  if (!model_ready)
    return Output{};

  //  assert(equalMatrix("Entering Solve()"));

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

  // comment the following line to use default starting point
  ComputeStartingPoint();

  // initialize residuals
  Residuals Res(m, n);
  ComputeResiduals_1234(Res);

  // initialize infeasibilities and mu
  double primal_infeas = Norm2(Res.res1) / Norm2(model.rhs);
  double dual_infeas = Norm2(Res.res4) / Norm2(model.obj);
  double mu = ComputeMu();
  double objective_value{};

  // ------------------------------------------
  // ---- MAIN LOOP ---------------------------
  // ------------------------------------------

  printf("\n");
  while (iter < options_.iteration_limit) {

    // Check that iterate is not NaN
    assert(!It.isNaN());

    // Stopping criterion
    if (iter > 0 && mu < options_.ipm_tolerance && // complementarity measure
        primal_infeas < options_.ipm_tolerance &&  // primal feasibility
        dual_infeas < options_.ipm_tolerance       // dual feasibility
    ) {
      printf("\n===== Optimal solution found =====\n\n");
      printf("Objective: %20.10e\n\n", DotProd(It.x, model.obj));
      break;
    }

    // Possibly print header
    if (iter % 20 == 0)
      printf(" iter            obj_v       pinf       dinf         mu        "
             "alpha_p    alpha_d\n");

    ++iter;

    // Compute diagonal scaling
    std::vector<double> scaling(n, 0.0);
    ComputeScaling(scaling);

    // Clear any existing INVERT now that scaling has changed
    invert.clear();

    //    printf("grep %d & %3.1f & %3.1f & %5.1g & %5.1g & %5.1g & %5.1g\n",
    //	   iter, It.x[0],It.x[1],1/scaling[0],1/scaling[1],1/scaling[2],1/scaling[3]);

    // Initialize Newton direction
    NewtonDir Delta(m, n);

    if (options_.predcor == 0) {

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
      SolveNewtonSystem(model.highs_a, scaling, Res, isCorrector, Delta);

      // Compute full Newton direction
      RecoverDirection(Res, isCorrector, Delta);

    } else {

      // *********************************************************************
      // PREDICTOR
      // *********************************************************************
      bool isCorrector = false;

      // Compute last two residuals for predictor
      ComputeResiduals_56(0, Delta, isCorrector, Res);

      // Solve Newton system for predictor
      SolveNewtonSystem(model.highs_a, scaling, Res, isCorrector, Delta);

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
      SolveNewtonSystem(model.highs_a, scaling, Res, isCorrector, DeltaCor);

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
    primal_infeas = Norm2(Res.res1) / Norm2(model.rhs);
    dual_infeas = Norm2(Res.res4) / Norm2(model.obj);
    objective_value = DotProd(It.x, model.obj);
    printf("%5d %20.10e %10.2e %10.2e %10.2e %10.2f %10.2f\n", iter,
           objective_value, primal_infeas, dual_infeas, mu, alpha_primal,
           alpha_dual);
  }

  // output struct
  Output out{};
  out.It = It;
  out.iterations = iter;
  out.primal_infeas = primal_infeas;
  out.dual_infeas = dual_infeas;
  out.mu = mu;

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

  // Check for NaN
  assert(!Res.isNaN());
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

void IPM_caller::reportOptions() {

  std::cout << "Option settings\n";
  std::cout << "NLA option:                  ";
  if (options_.nla == kOptionNlaCg) {
    std::cout << "Newton CG\n";
  } else if (options_.nla == kOptionNlaAugmented || options_.nla == kOptionNlaAugmentedCg) {
    std::cout << "Augmented direct\n";
  } else if (options_.nla == kOptionNlaNewton || options_.nla == kOptionNlaNewtonCg) {
    std::cout << "Newton direct\n";
  }
  std::cout << "Predictor-corrector:         ";
  if (this->options_.predcor) {
    std::cout << "On\n";
  } else {
    std::cout << "Off\n";
  }
  std::cout << "Max number of dense columns: " << this->options_.max_dense_col << "\n";
  std::cout << "Tolerance for dense columns: " << this->options_.dense_col_tolerance << "\n";
  std::cout << "IPM iteration limit:         " << this->options_.iteration_limit << "\n";
  std::cout << "IPM tolerance:               " << this->options_.ipm_tolerance << "\n";
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

  // Identify whether CG should be used
  const bool use_cg = options_.nla == kOptionNlaCg ||
                      options_.nla == kOptionNlaAugmentedCg ||
                      options_.nla == kOptionNlaNewtonCg;
  // Identify whether the augmented system should be solved directly
  const bool use_direct_augmented =
      options_.nla == kOptionNlaAugmented || options_.nla == kOptionNlaAugmentedCg;
  // Identify whether the Newton system should be solved directly
  const bool use_direct_newton =
      options_.nla == kOptionNlaNewton || options_.nla == kOptionNlaNewtonCg;
  // Shouldn't be trying to solve both the augmented and Newton system!
  assert(!use_direct_augmented || !use_direct_newton);

  // Until the direct solvers are implemented correctly, the IPM
  // solver is driven by CG, so check for this
  //  assert(use_cg);

  // Identify whether the CG result should be used to check the result
  // obtained using the direct solver
  const bool check_with_cg = use_cg && options_.nla != kOptionNlaCg;

  // Augmented system is
  //
  // [-scaling A^T][dx] = [res7]
  // [A        0  ][dy]   [res1]
  //
  std::vector<double> theta;
  scaling2theta(scaling, theta);

  const bool first_call_with_theta = !invert.valid;
  ExperimentData experiment_data;
  experiment_data.reset();
  if (first_call_with_theta) experiment_data.analyseTheta(theta, false);
  if (use_cg || use_direct_newton) {
    // Have to solve the Newton system
    //
    // A.scaling.A^T Delta.y = res8
    //

    // Compute res8
    std::vector<double> res8{
        ComputeResiduals_8(highs_a, scaling, Res, res7, isCorrector)};

    // Solve normal equations
    //
    // Currently this is done using Conjugate Gradient. The solution for
    // Delta.y can be substituted with a positive definite factorization.
    //
    if (use_cg) {
      NormalEquations N(highs_a, scaling);
      CG_solve(N, res8, kCgTolerance, kCgIterationLimit, Delta.y, nullptr);
      const bool compute_residual_error = false;
      if (compute_residual_error) {
        std::pair<double, double> residual_error =
            residualErrorNewton(highs_a, theta, res8, Delta.y);
        if (residual_error.first > 1e-12)
          printf("CG solve abs (rel) residual error = %g (%g)\n",
                 residual_error.first, residual_error.second);
      }
    }
    if (use_direct_newton) {
      // Solve the Newton system directly into newton_delta_y
      std::vector<double> newton_delta_y;
      newton_delta_y.assign(m, 0);
      if (first_call_with_theta) {
        int newton_invert_status =
            newtonInvert(highs_a, theta, invert, options_.max_dense_col,
                         options_.dense_col_tolerance, experiment_data);
        if (newton_invert_status) return newton_invert_status;
      } else {
	// Just set this to avoid assert when writing out
	// experiment_data in the event of a solution error
	experiment_data.system_type = kSystemTypeNewton;
      }
      int newton_solve_status = newtonSolve(
          highs_a, theta, res8, newton_delta_y, invert, experiment_data);
      if (first_call_with_theta) 
        experiment_data.condition = newtonCondition(highs_a, theta, invert);
      if (newton_solve_status) return newton_solve_status;
      if (check_with_cg) {
        double inf_norm_solution_diff = infNormDiff(newton_delta_y, Delta.y);
        if (inf_norm_solution_diff > kSolutionDiffTolerance) {
          std::cout << "Newton Direct-CG solution error = "
                    << inf_norm_solution_diff << "\n";
          std::cout << experiment_data << "\n";
          assert(111 == 333);
        }
      }
      if (!use_cg) {
        // Once CG is not being used, use the Newton solution for dy
        Delta.y = newton_delta_y;
      }
    }

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
    // Solve augmented system directly
    std::vector<double> augmented_delta_x;
    augmented_delta_x.assign(n, 0);
    std::vector<double> augmented_delta_y;
    augmented_delta_y.assign(m, 0);
    if (first_call_with_theta) {
      int augmented_invert_status =
          augmentedInvert(highs_a, theta, invert, experiment_data);
      if (augmented_invert_status) return augmented_invert_status;
    } else {
      // Just set this to avoid assert when writing out
      // experiment_data in the event of a solution error
      experiment_data.system_type = kSystemTypeNewton;
    }
    // When solving the augmented system, the right side for the
    // predictor should be (res7; res1), but for the corrector it
    // should be (res7; 0)
    if (isCorrector) {
      // Check that res1 has been zeroed
      const double res1_norm = Norm2(Res.res1);
      assert(!res1_norm);
    }
    augmentedSolve(highs_a, theta, res7, Res.res1, augmented_delta_x,
                   augmented_delta_y, invert, experiment_data);
    experiment_data.nla_time.total += experiment_data.nla_time.solve;
    if (first_call_with_theta) 
      experiment_data.condition = augmentedCondition(highs_a, theta, invert);
    if (check_with_cg) {
      double inf_norm_solution_diff = infNormDiff(augmented_delta_x, Delta.x);
      if (inf_norm_solution_diff > kSolutionDiffTolerance) {
        std::cout << "Augmented Direct-CG solution error for Delta.x = "
                  << inf_norm_solution_diff << "\n";
        std::cout << experiment_data << "\n";
        assert(111 == 333);
      }
      inf_norm_solution_diff = infNormDiff(augmented_delta_y, Delta.y);
      if (inf_norm_solution_diff > kSolutionDiffTolerance) {
        std::cout << "Augmented Direct-CG solution error for Delta.y = "
                  << inf_norm_solution_diff << "\n";
        std::cout << experiment_data << "\n";
        assert(111 == 333);
      }
    }
    if (!use_cg) {
      // Once CG is not being used, use the augmented solution for dx and dy
      Delta.x = augmented_delta_x;
      Delta.y = augmented_delta_y;
    }
  }
  if (first_call_with_theta)
    experiment_data_record.push_back(experiment_data);
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

  // Check for NaN
  assert(!Delta.isNaN());
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

  // use y to store b-A*x
  It.y = model.rhs;
  model.highs_a.alphaProductPlusY(-1.0, It.x, It.y);

  // solve A*A^T * dx = b-A*x with CG and store the result in temp_m
  std::vector<double> temp_scaling(n, 1.0);
  NormalEquations N(model.highs_a, temp_scaling);

  std::vector<double> temp_m(m);
  int cg_iter{};
  CG_solve(N, It.y, 1e-4, 100, temp_m, &cg_iter);
  int cg_iter_starting_point{cg_iter};

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
    It.xl[i] = It.x[i] - model.lower[i];
    It.xu[i] = model.upper[i] - It.x[i];

    violation = std::min(violation, It.xl[i]);
    violation = std::min(violation, It.xu[i]);
  }

  // shift to be positive
  violation = 1.0 + std::max(0.0, -1.5 * violation);
  VectorAdd(It.xl, violation);
  VectorAdd(It.xu, violation);
  // *********************************************************************

  // *********************************************************************
  // y starting point
  // *********************************************************************
  // compute A*c
  std::fill(temp_m.begin(), temp_m.end(), 0.0);
  model.highs_a.alphaProductPlusY(1.0, model.obj, temp_m);

  // compute (A*A^T)^{-1} * A*c and store in y
  CG_solve(N, temp_m, 1e-4, 100, It.y, &cg_iter);
  cg_iter_starting_point += cg_iter;
  printf("Starting point required %d CG iterations\n", cg_iter_starting_point);
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
