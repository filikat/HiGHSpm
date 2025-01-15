#ifndef IPM_H
#define IPM_H

#include <string>

#include "../FactorHiGHS/FactorHiGHS.h"
#include "CgSolver.h"
#include "FactorHiGHSSolver.h"
#include "IpmModel.h"
#include "Ipm_aux.h"
#include "Ipm_const.h"
#include "LinearSolver.h"
#include "VectorOperations.h"
#include "util/HighsSparseMatrix.h"

class Ipm {
  // LP model
  IpmModel model_;

  // Objects used during iterations
  Iterate it_{};
  Residuals res_{};
  NewtonDir delta_{};

  // Linear solver interface
  std::unique_ptr<LinearSolver> LS_;

  // Size of the problem
  int m_{};
  int n_{};

  // Iterations counters
  int iter_{};
  int bad_iter_{};

  // Indicators
  double primal_infeas_{};
  double dual_infeas_{};
  double primal_obj_{};
  double dual_obj_{};
  double pd_gap_{};
  double mu_{};

  // Other statistics
  double min_prod_{};
  double max_prod_{};

  // Theta^-1
  std::vector<double> scaling_{};

  // Stepsizes
  double alpha_primal_{};
  double alpha_dual_{};

  // Coefficient for reduction of mu
  double sigma_{};

  // Status of the solver
  std::string ipm_status_ = "Max iter";

  // Run-time options
  Options options_{};

  // Timer for iterations
  Clock clock_;

 public:
  // ===================================================================================
  // Load an LP:
  //
  //  min   obj^T * x
  //  s.t.  Ax {<=,=,>=} rhs
  //        lower <= x <= upper
  //
  // Transform constraints in equalities by adding slacks to inequalities:
  //  <= : add slack    0 <= s_i <= +inf
  //  >= : add slack -inf <= s_i <=    0
  // ===================================================================================
  void load(const int num_var,           // number of variables
            const int num_con,           // number of constraints
            const double* obj,           // objective function c
            const double* rhs,           // rhs vector b
            const double* lower,         // lower bound vector
            const double* upper,         // upper bound vector
            const int* A_ptr,            // column pointers of A
            const int* A_rows,           // row indices of A
            const double* A_vals,        // values of A
            const int* constraints,      // type of constraints
            const std::string& pb_name,  // problem name
            const Options& options       // options
  );

  // ===================================================================================
  // Solve the LP
  // ===================================================================================
  Output solve();

 private:
  // ===================================================================================
  // Compute:
  //
  //  mu = \sum xl(i) * zl(i) + \sum xu(j) * zu(j)
  //
  // for variables: i with a finite lower bound
  //                j with a finite upper bound
  // ===================================================================================
  void computeMu();

  // ===================================================================================
  // Compute:
  //
  //  res1 = rhs - A * x
  //  res2 = lower - x + xl
  //  res3 = upper - x - xu
  //  res4 = c - A^T * y - zl + zu
  //
  // Components of residuals 2,3 are set to zero if the corresponding
  // upper/lower bound is not finite.
  // ===================================================================================
  void computeResiduals1234();

  // ===================================================================================
  // Compute:
  //
  //  res5 = sigma * mu * e - Xl * Zl * e
  //  res6 = sigma * mu * e - Xu * Zu * e
  //
  //  Or
  //
  //  res5 = sigma * mu * e - D_aff Xl * D_aff Zl * e
  //  res6 = sigma * mu * e - D_aff Xu * D_aff Zu * e
  //
  // Components of residuals 5,6 are set to zero if the corresponding
  // upper/lower bound is not finite.
  // ===================================================================================
  void computeResiduals56();

  // ===================================================================================
  // Compute:
  //
  //  res7 = res4 - Xl^{-1} * (res5 + Zl * res2) + Xu^{-1} * (res6 - Zu * res3)
  //
  // (the computation of res7 takes into account only the components for which
  // the correspoding upper/lower bounds are finite)
  // ===================================================================================
  std::vector<double> computeResiduals7();

  // ===================================================================================
  // Compute:
  //
  //  res8 = res1 + A * Theta * res7
  // ===================================================================================
  std::vector<double> computeResiduals8(const std::vector<double>& res7);

  // ===================================================================================
  // Compute diagonal scaling Theta^{-1}
  //
  //  Theta^{-1}_{ii} = zl(i) / xl(i) + zu(i) / xu(i)
  //
  // Theta^{-1} only considers the terms above if the corresponding upper/lower
  // bound is finite.
  // ===================================================================================
  void computeScaling();

  // ===================================================================================
  // Solve:
  //
  // ___Augmented system___
  //
  //      [ -Theta^{-1}  A^T ] [ Deltax ] = [ res7 ]
  //      [ A            0   ] [ Deltay ] = [ res1 ]
  //
  // with:
  //  res7 = res4 - Xl^{-1} * (res5 + Zl * res2) + Xu^{-1} * (res6 - Zu * res3)
  //  Theta^{-1} = diag( scaling )
  //
  // (the computation of res7 takes into account only the components for which
  // the correspoding upper/lower bounds are finite)
  //
  // OR
  //
  // ___Normal equations___
  //
  //      A * Theta * A^T * Deltay = res8
  //      Delta x = Theta * (A^T* Deltay - res7)
  //
  // with:
  //  res8 = res1 + A * Theta * res7
  // ===================================================================================
  bool solveNewtonSystem(NewtonDir& delta);

  // ===================================================================================
  // Reconstruct the solution of the full Newton system:
  //
  //  Deltaxl = Deltax - res2
  //  Deltaxu = res3 - Deltax
  //  Deltazl = Xl^{-1} * (res5 - zl * Deltaxl)
  //  Deltazu = Xu^{-1} * (res6 - zu * Deltaxu)
  // ===================================================================================
  bool recoverDirection(NewtonDir& delta);

  // ===================================================================================
  // Find alpha_primal:
  //
  //  x  + alpha_primal * Deltax
  //  xl + alpha_primal * Deltaxl > 0     (if lower bound finite)
  //  xu + alpha_primal * Deltaxu > 0     (if upper bound finite)
  //
  // Find alpha_dual:
  //
  //  y  + alpha_dual * Deltay
  //  zl + alpha_dual * Deltazl > 0       (if lower bound finite)
  //  zu + alpha_dual * Deltazu > 0       (if upper bound finite)
  //
  // Step-sizes are scaled down by kInteriorScaling < 1, to guarantee that no
  // component of the new iterate is equal to zero.
  //
  // If corrector is valid, the direction used is delta + weight * corrector
  // ===================================================================================
  void computeStepSizes(double& alpha_primal, double& alpha_dual,
                        const NewtonDir& delta,
                        const NewtonDir* corrector = nullptr,
                        double weight = 1.0) const;

  // ===================================================================================
  // Make the step in the Newton direction with appropriate stepsizes.
  // ===================================================================================
  void makeStep();

  // ===================================================================================
  // Compute the Mehrotra starting point.
  // In doing so, CG is used to solve two linear systems with matrix A*A^T.
  // This task does not need a factorization and can continue to use CG, because
  // these linear systems are very easy to solve with CG.
  // ===================================================================================
  void computeStartingPoint();

  // ===================================================================================
  // Given the predictor direction, compute predicted mu
  //    mu_aff = (xl + alpha_p * DeltaAff xl)' * (zl + alpha_d * DeltaAff zl) +
  //             (xu + alpha_p * DeltaAff xu)' * (zu + alpha_d * DeltaAff zu)
  //    mu_aff /= num_finite_bounds
  //
  // and return sigma for the corrector direction
  //    sigma = ( mu_aff / mu )^3
  // ===================================================================================
  void computeSigma();

  // ===================================================================================
  // Compute the residuals for the computation of multiple centrality
  // correctors.
  // ===================================================================================
  void computeResidualsMcc();

  // ===================================================================================
  // Iteratively compute correctors, until they improve the stepsizes.
  // Based on Gondzio, "Multiple centrality corrections in a primal-dual method
  // for linear programming" and Colombo, Gondzio, "Further Development of
  // Multiple Centrality Correctors for Interior Point Methods".
  // ===================================================================================
  bool centralityCorrectors();

  // ===================================================================================
  // Given the current direction delta and the latest corrector, compute the
  // best primal and dual weights, that maximize the primal and dual stepsize.
  // ===================================================================================
  void computeBestWeight(const NewtonDir& delta, const NewtonDir& corrector,
                         double& wp, double& wd, double& alpha_p,
                         double& alpha_d) const;

  // ===================================================================================
  // Compute
  // - primal infeasibility
  // - dual infeasibility
  // - primal objective
  // - dual objective
  // - primal-dual relative gap
  // ===================================================================================
  void computeIndicators();

  // ===================================================================================
  // Compute the complementarity products and potentially alter them.
  // ===================================================================================
  void computeProducts();

  // ===================================================================================
  // If the current iterate is nan or inf, abort the iterations.
  // ===================================================================================
  bool checkIterate();

  // ===================================================================================
  // If too many bad iterations happened consecutively, abort the iterations.
  // ===================================================================================
  bool checkBadIter();

  // ===================================================================================
  // Check the termination criterion:
  //  - primal infeasibility < tolerance
  //  - dual infeasiblity    < tolerance
  //  - relative dual gap    < tolerance
  // ===================================================================================
  bool checkTermination();

  void printInfo() const;
  void printHeader() const;
  void printOutput() const;
  void collectData() const;
};

#endif
