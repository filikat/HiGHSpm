#ifndef IPM_CALLER
#define IPM_CALLER

#include "IPM_model.h"
#include "IPM_aux.h"
#include "SparseMatrix.h"
#include "VectorOperations.h"
#include "NormalEquations.h"
#include "ConjugateGradient.h"
#include "IPM_aux.h"



class IPM_caller{

    IPM_model model;
    int m{};
    int n{};
    bool model_ready{false};

    // IPM parameters
    const double sigma_i = 0.5;
    const double sigma_min = 0.05;
    const double sigma_max = 0.95;
    const double alpha_interior_scaling = 0.99;

    // IPM iterate
    Iterate It{};



public:


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
    void Load(      // - - - - - - - - - - - - - - -  length
        //INPUT
        const int num_var,      // number of variables
        const int num_con,      // number of constraints
        const double* obj,      // objective function c           num_var
        const double* rhs,      // rhs vector b                   num_con
        const double* lower,    // lower bound vector             num_var
        const double* upper,    // upper bound vector             num_var
        const int* A_colptr,    // column pointers of A           num_var + 1
        const int* A_rowind,    // row indices of A               A_colptr[num_var]
        const double* A_values, // values of A                    A_colptr[num_var]
        const int* constraints  // type of constraints            num_con
                                // 1: >=, 0: =, -1: <=
    );


// ===================================================================================
// SOLVE THE LP
// ===================================================================================
    Output Solve(
        //INPUT
        const int itmax,        // maximum number of IPM iterations
        const double tol        // tolerance on primal/dual feasibility and mu
    );


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
// COMPUTE RESIDUALS 
// ===================================================================================
// Compute:
//
//  res1 = rhs - A * x
//  res2 = lower - x + xl
//  res3 = upper - x - xu
//  res4 = c - A^T * y - zl + zu
//  res5 = sigma * mu * e - Xl * Zl * e
//  res6 = sigma * mu * e - Xu * Zu * e
//
// Components of residuals 2,3,5,6 are set to zero if the corresponding upper/lower
// bound is not finite.
//
// ===================================================================================
    void ComputeResiduals(
        //INPUT
        const double sigmaMu,   // sigma * mu
        //OUTPUT
        Residuals& Res          // residuals
    );


// ===================================================================================
// COMPUTE SCALING
// ===================================================================================
// Compute diagonal scaling Theta^{-1}
//
//  Theta^{-1}_{ii} = zl(i) / xl(i) + zu(i) / xu(i)
//
// Theta^{-1} only considers the terms above if the corresponding upper/lower
// bound are finite.
//
// ===================================================================================
    void ComputeScaling(
        //OUTPUT
        std::vector<double>& scaling        // diagonal scaling, length n
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
// (the computation of res7 takes into account only the components for which the
//  correspoding upper/lower bounds are finite)
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
    void SolveNewtonSystem(
        // INPUT
        const SparseMatrix& A,                  // constraint matrix
        const std::vector<double>& scaling,     // diagonal scaling, length n
        const Residuals& Res,                   // current residuals
        // OUTPUT
        NewtonDir& Delta                        // Newton direction
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
        //INPUT
        const Residuals& Res,                   // current residuals
        //OUTPUT
        NewtonDir& Delta                        // Newton directions
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
// Step-sizes are scaled down by alpha_interior_scaling < 1, to guarantee that no
// component of the new iterate is equal to zero.
//
// ===================================================================================                  
    void ComputeStepSizes(
        //INPUT            
        const NewtonDir& Delta,                 // Newton direction
        //OUTPUT
        double& alpha_primal,                   // primal step-size
        double& alpha_dual                      // dual step-size
    );

};

#endif