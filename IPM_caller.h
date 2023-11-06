#ifndef IPM_CALLER
#define IPM_CALLER

#include "ConjugateGradient.h"
#include "Direct.h"
#include "IPM_aux.h"
#include "IPM_const.h"
#include "IPM_model.h"
#include "Metis_caller.h"
#include "NormalEquations.h"
#include "VectorOperations.h"
#include "util/HighsSparseMatrix.h"

class IPM_caller {
  IPM_model model;

  int m{};
  int n{};
  bool model_ready{false};
  IpmInvert invert;

  // IPM parameters
  const double sigma_i = 0.5;
  const double sigma_min = 0.05;
  const double sigma_max = 0.95;
  const double alpha_interior_scaling = 0.99;

  // IPM iterate
  Iterate It{};

  Metis_caller Metis_data;

 public:
  // ===================================================================================
  // Run-time options
  // ===================================================================================
  int option_nla = kOptionNlaDefault;
  int option_max_dense_col = kOptionMaxDenseColDefault;
  double option_dense_col_tolerance = kOptionDenseColToleranceDefault;
  int option_predcor = kOptionPredCorDefault;

  // Direct solver experiment data record
  std::vector<ExperimentData> experiment_data_record;

  // ===================================================================================
  // LOAD THE PROBLEM
  // ===================================================================================
  // Load as input an LP:
  //
  //  min   obj^T * x
  //  s.t.  Ax {<=,=,>=} rhs
  //        lower <= x <= upper
  //
  // Transform constraints in equalities by adding slacks to inequalities:
  //  <= : add slack    0 <= s_i <= +inf
  //  >= : add slack -inf <= s_i <=    0
  //
  // ===================================================================================
  void Load(                // - - - - - - - - - - - - - - -  length
                            // INPUT
      const int num_var,    // number of variables
      const int num_con,    // number of constraints
      const double* obj,    // objective function c           num_var
      const double* rhs,    // rhs vector b                   num_con
      const double* lower,  // lower bound vector             num_var
      const double* upper,  // upper bound vector             num_var
      const int* A_colptr,  // column pointers of A           num_var + 1
      const int* A_rowind,  // row indices of A               A_colptr[num_var]
      const double* A_values,  // values of A A_colptr[num_var]
      const int* constraints   // type of constraints            num_con
                               // 1: >=, 0: =, -1: <=
  );

  // ===================================================================================
  // SOLVE THE LP
  // ===================================================================================
  Output Solve();

 private:
  // ===================================================================================
  // COMPUTE MU
  // ===================================================================================
  // Compute:
  //
  //  mu = \sum xl(i) * zl(i) + \sum xu(j) * zu(j)
  //
  // for variables: i with a finite lower bound
  //                j with a finite upper bound
  //
  // ===================================================================================
  double ComputeMu();

  // ===================================================================================
  // COMPUTE RESIDUALS 1,2,3,4
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
  //
  // ===================================================================================
  void ComputeResiduals_1234(
      // OUTPUT
      Residuals& Res  // residuals
  );

  // ===================================================================================
  // COMPUTE RESIDUALS 5,6
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
  //
  // ===================================================================================
  void ComputeResiduals_56(
      // INPUT
      const double sigmaMu,       // sigma * mu
      const NewtonDir& DeltaAff,  // affine scaling direction (if corrector)
      bool isCorrector,           // true if corrector, false if predictor
      // OUTPUT
      Residuals& Res  // residuals
  );

  // ===================================================================================
  // COMPUTE RESIDUAL 7
  // ===================================================================================
  // Compute:
  //
  //  res7 = res4 - Xl^{-1} * (res5 + Zl * res2) + Xu^{-1} * (res6 - Zu * res3)
  //
  //  Or
  //
  //  res7 = - Xl^{-1} * res5 + Xu^{-1} * res6
  //
  // (the computation of res7 takes into account only the components for which
  // the correspoding upper/lower bounds are finite)
  //
  // ===================================================================================
  std::vector<double> ComputeResiduals_7(
      // INPUT
      const Residuals& Res,     // residuals
      bool isCorrector = false  // true if corrector, false is predictor
  );

  // ===================================================================================
  // COMPUTE RESIDUAL 8
  // ===================================================================================
  // Compute:
  //
  //  res8 = res1 + A * Theta * res7
  //
  //  Or
  //
  //  res8 = A * Theta * res7
  //
  // ===================================================================================
  std::vector<double> ComputeResiduals_8(
      // INPUT
      const HighsSparseMatrix& highs_a,    // constraint matrix
      const std::vector<double>& scaling,  // scaling vector
      const Residuals& Res,                // residuals
      const std::vector<double>& res7,     // residual 7
      bool isCorrector = false  // true if corrector, false is predictor
  );

  // ===================================================================================
  // COMPUTE SCALING
  // ===================================================================================
  // Compute diagonal scaling Theta^{-1}
  //
  //  Theta^{-1}_{ii} = zl(i) / xl(i) + zu(i) / xu(i)
  //
  // Theta^{-1} only considers the terms above if the corresponding upper/lower
  // bound is finite.
  //
  // ===================================================================================
  void ComputeScaling(
      // OUTPUT
      std::vector<double>& scaling  // diagonal scaling, length n
  );

  // ===================================================================================
  // SOLVE NEWTON SYSTEM
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
  //
  //
  // ___WARNING___
  //  scaling may contain zeros, if there are free variables.
  //  Do not form Theta = diag( scaling )^{-1} in this case.
  //
  // ===================================================================================
  int SolveNewtonSystem(
      // INPUT
      const HighsSparseMatrix& highs_a,    // constraint matrix
      const std::vector<double>& scaling,  // diagonal scaling, length n
      const Residuals& Res,                // current residuals
      bool isCorrector,  // true if corrector, false if predictor
      // OUTPUT
      NewtonDir& Delta  // Newton direction
  );

  // ===================================================================================
  // FULL NEWTON DIRECTION
  // ===================================================================================
  // Reconstruct the solution of the full Newton system:
  //
  //  Deltaxl = Deltax - res2
  //  Deltaxu = res3 - Deltax
  //  Deltazl = Xl^{-1} * (res5 - zl * Deltaxl)
  //  Deltazu = Xu^{-1} * (res6 - zu * Deltaxu)
  //
  // ===================================================================================
  void RecoverDirection(
      // INPUT
      const Residuals& Res,  // current residuals
      bool isCorrector,      // true if corrector, false if predictor
      // OUTPUT
      NewtonDir& Delta  // Newton directions
  );

  // ===================================================================================
  // COMPUTE STEP-SIZES
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
  // Step-sizes are scaled down by alpha_interior_scaling < 1, to guarantee that
  // no component of the new iterate is equal to zero.
  //
  // ===================================================================================
  void ComputeStepSizes(
      // INPUT
      const NewtonDir& Delta,  // Newton direction
      // OUTPUT
      double& alpha_primal,  // primal step-size
      double& alpha_dual     // dual step-size
  );

  // ===================================================================================
  // COMPUTE STARTING POINT
  // ===================================================================================
  // Compute the Mehrotra starting point.
  // In doing so, CG is used to solve two linear systems with matrix A*A^T.
  // This task does not need a factorization and can continue to use CG, because
  // these linear systems are very easy to solve with CG.
  //
  // ===================================================================================
  void ComputeStartingPoint();

  // ===================================================================================
  // COMPUTE SIGMA FOR CORRECTOR
  // ===================================================================================
  //  Given the predictor direction, compute predicted mu
  //    mu_aff = (xl + alpha_p * DeltaAff xl)' * (zl + alpha_d * DeltaAff zl) +
  //             (xu + alpha_p * DeltaAff xu)' * (zu + alpha_d * DeltaAff zu)
  //    mu_aff /= num_finite_bounds
  //
  //  and return sigma for the corrector direction
  //    sigma = ( mu_aff / mu )^3
  //
  // ===================================================================================
  double ComputeSigmaCorrector(
      // INPUT
      const NewtonDir& DeltaAff,  // Predictor Newton direction
      double mu                   // mu of previous iteration
  );

  void CheckResiduals(const NewtonDir& Delta, const Residuals& Res) const;
};

#endif
