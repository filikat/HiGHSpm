#include <iostream>
#include <cmath>
#include "IPM_caller.h"

#include <fstream>


// =======================================================================
// LOAD THE PROBLEM
// =======================================================================
void IPM_caller::Load(
    const int num_var,
    const int num_con,
    const double* obj,
    const double* rhs,
    const double* lower,
    const double* upper,
    const int* A_colptr,
    const int* A_rowind,
    const double* A_values,
    const int* constraints
){

    if (!obj || !rhs || !lower || !upper || !A_colptr || 
        !A_rowind || !A_values || !constraints) return;

    // count how many slacks are needed
    int num_slacks{};
    for (int i=0; i<num_con; ++i){
        if (constraints[i] != 0){

            ++num_slacks;

            if (constraints[i] != 1 && constraints[i] != -1){
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

    for (int i=0; i<num_var; ++i){

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

    for (int i=0; i<A_nnz; ++i){

        // copy original row indices
        temp_rowind[i] = A_rowind[i];

        // copy original values
        temp_values[i] = A_values[i];

    }

    // index of the current slack
    int slack_ind{num_var};

    for (int i=0; i<num_con; ++i){

        // copy vector b
        model.rhs[i] = rhs[i];

        // if constraint is inequality, add a slack variable
        if (constraints[i] != 0){

            // lower/upper bound for new slack
            if (constraints[i] == 1){
                model.lower[slack_ind] = - INF;
                model.upper[slack_ind] = 0.0;
            }else{
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

    model.A = SparseMatrix(temp_rowind,temp_colptr,temp_values,num_con,num_var + num_slacks);

    m = model.num_con;
    n = model.num_var;
    model_ready = true;

}


// =======================================================================
// SOLVE THE LP
// =======================================================================
Output IPM_caller::Solve(const int itmax, const double tol){

    // solve only if model is loaded
    if (!model_ready) return Output{};

    // ------------------------------------------
    // ---- INITIALIZE --------------------------
    // ------------------------------------------

    // iterations counter
    int iter{};

    // initialize stepsize
    double alpha_primal{};
    double alpha_dual{};

    // initialize starting point
    It = Iterate(m,n);

    // initialize residuals
    Residuals Res(m,n);

    // initialize infeasibilities and mu
    double primal_infeas{};
    double dual_infeas{};
    double mu{};

    // print header
    printf("\n\n  iter     pinf      dinf        mu        alpha_p    alpha_d\n");


    // ------------------------------------------
    // ---- MAIN LOOP ---------------------------
    // ------------------------------------------
    while (iter < itmax){

        // Stopping criterion
        if (iter > 0 &&                                 
            mu < tol &&                                 // complementarity measure
            Norm2(Res.res1)/Norm2(model.rhs) < tol &&   // primal feasibility
            Norm2(Res.res4)/Norm2(model.obj) < tol      // dual feasibility
            ){
                printf("\n===== Optimal solution found =====\n\n");
                printf("Objective: %20.10e\n\n",DotProd(It.x,model.obj));
                break;
            }

        ++iter;

        // Heuristic to choose sigma
        // (this is redundant once we implement predictor-corrector)
        double sigma{};
        if (iter == 1){
            sigma = sigma_i;
        }else{
            sigma = pow(std::max(1.0-alpha_primal,1.0-alpha_dual), 5.0);
            sigma = std::max(sigma,sigma_min);
            sigma = std::min(sigma,sigma_max);
        }

        // Find mu
        mu = ComputeMu();

        // Compute residuals
        ComputeResiduals(sigma*mu,Res);

        // Compute diagonal scaling
        std::vector<double> scaling(n,0.0);
        ComputeScaling(scaling);

        // Initialize Newton direction
        NewtonDir Delta(m,n);

        // Solve Newton system
        SolveNewtonSystem(model.A,scaling,Res,Delta);

        // Compute full Newton direction        
        RecoverDirection(Res,Delta);

        // Find step-sizes
        ComputeStepSizes(Delta,alpha_primal,alpha_dual);

        // Make the step
        VectorAdd(It.x,Delta.x,alpha_primal);
        VectorAdd(It.xl,Delta.xl,alpha_primal);
        VectorAdd(It.xu,Delta.xu,alpha_primal);
        VectorAdd(It.y,Delta.y,alpha_dual);
        VectorAdd(It.zl,Delta.zl,alpha_dual);
        VectorAdd(It.zu,Delta.zu,alpha_dual);

        // Print output to screen
        primal_infeas = Norm2(Res.res1) / Norm2(model.rhs);
        dual_infeas = Norm2(Res.res4) / Norm2(model.obj);
        printf("%5d %10.2e %10.2e %10.2e %10.2f %10.2f\n",iter,primal_infeas,dual_infeas,mu,alpha_primal,alpha_dual);

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
double IPM_caller::ComputeMu(){
    double mu = 0.0;
    int number_finite_bounds{};
    for (int i=0; i<n; ++i){
        if (model.has_lb(i)){
            mu += It.xl[i] * It.zl[i];
            ++number_finite_bounds;
        }
        if (model.has_ub(i)){
            mu += It.xu[i] * It.zu[i];
            ++number_finite_bounds;
        }
    }
    return (mu / number_finite_bounds);
}


// =======================================================================
// COMPUTE RESIDUALS 
// =======================================================================
void IPM_caller::ComputeResiduals(const double sigmaMu,Residuals& Res){
    
    // res1
    Res.res1 = model.rhs;
    mat_vec(model.A,It.x,Res.res1,-1.0,'n');

    // res2
    for (int i=0; i<n; ++i){
        if (model.has_lb(i)){
            Res.res2[i] = model.lower[i] - It.x[i] + It.xl[i];
        }else{
            Res.res2[i] = 0.0;
        }
    }

    // res3
    for (int i=0; i<n; ++i){
        if (model.has_ub(i)){
            Res.res3[i] = model.upper[i] - It.x[i] - It.xu[i];
        }else{
            Res.res3[i] = 0.0;
        }
    }

    // res4
    Res.res4 = model.obj;
    mat_vec(model.A,It.y,Res.res4,-1.0,'t');
    for (int i=0; i<n; ++i){
        if (model.has_lb(i)){
            Res.res4[i] -= It.zl[i];
        }
        if (model.has_ub(i)){
            Res.res4[i] += It.zu[i];
        }
    }

    // res5
    for (int i=0; i<n; ++i){
        if (model.has_lb(i)){
            Res.res5[i] = sigmaMu - It.xl[i] * It.zl[i];
        }else{
            Res.res5[i] = 0.0;
        }
    }

    // res6
    for (int i=0; i<n; ++i){
        if (model.has_ub(i)){
            Res.res6[i] = sigmaMu - It.xu[i] * It.zu[i];
        }else{
            Res.res6[i] = 0.0;
        }
    }

}


// =======================================================================
// COMPUTE SCALING
// =======================================================================
void IPM_caller::ComputeScaling(std::vector<double>& scaling){

    for (int i=0; i<n; ++i){
        if (model.has_lb(i)){
            scaling[i] += It.zl[i] / It.xl[i];
        }
        if (model.has_ub(i)){
            scaling[i] += It.zu[i] / It.xu[i];
        }
    }
    
}


// =======================================================================
// SOLVE NEWTON SYSTEM
// =======================================================================
void IPM_caller::SolveNewtonSystem(
    const SparseMatrix& A,
    const std::vector<double>& scaling,
    const Residuals& Res,
    NewtonDir& Delta 
){

    // Compute res7
    // *********************************************************************
    std::vector<double> res7(Res.res4);

    for (int i=0; i<n; ++i){
        if (model.has_lb(i)){
            res7[i] -= ( (Res.res5[i] + It.zl[i] * Res.res2[i]) / It.xl[i] );
        }
        if (model.has_ub(i)){
            res7[i] += ( (Res.res6[i] - It.zu[i] * Res.res3[i]) / It.xu[i] );
        }
    }
    // *********************************************************************


    // Compute res8
    // *********************************************************************
    std::vector<double> res8(Res.res1);
    std::vector<double> temp(res7);

    // temp = Theta * res7
    VectorDivide(temp,scaling);

    // res8 += A * temp
    mat_vec(A,temp,res8,1.0,'n');
    temp.clear();
    // *********************************************************************


    // Solve normal equations
    // Currently this is done using Conjugate Gradient. The solution for 
    // Delta.y can be substituted with a positive definite factorization.
    NormalEquations N(A,scaling);
    CG_solve(N,res8,1e-12,5000,Delta.y);
    

    // Compute Delta.x
    // *********************************************************************
    // Deltax = A^T * Deltay - res7;
    Delta.x = res7;
    mat_vec(A,Delta.y,Delta.x,-1.0,'t');
    VectorScale(Delta.x,-1.0);

    // Deltax = Theta * Deltax
    VectorDivide(Delta.x,scaling);
    // *********************************************************************

}


// =======================================================================
// FULL NEWTON DIRECTION 
// =======================================================================
void IPM_caller::RecoverDirection(const Residuals& Res,NewtonDir& Delta){
    
    // Deltaxl
    Delta.xl = Delta.x;
    VectorAdd(Delta.xl,Res.res2,-1.0);

    // Deltaxu
    Delta.xu = Res.res3;
    VectorAdd(Delta.xu,Delta.x,-1.0);

    // Deltazl
    Delta.zl = Res.res5;
    VectorAddMult(Delta.zl,It.zl,Delta.xl,-1.0);
    VectorDivide(Delta.zl,It.xl);

    // Deltazu
    Delta.zu = Res.res6;
    VectorAddMult(Delta.zu,It.zu,Delta.xu,-1.0);
    VectorDivide(Delta.zu,It.xu);

}

// =======================================================================
// COMPUTE STEP-SIZES
// =======================================================================                  
void IPM_caller::ComputeStepSizes(
    const NewtonDir& Delta,
    double& alpha_primal,
    double& alpha_dual 
){

    alpha_primal = 1.0;
    for (int i=0; i<n; ++i){
        if (Delta.xl[i] < 0 && model.has_lb(i)){
            alpha_primal = std::min(alpha_primal, - It.xl[i] / Delta.xl[i]);
        }
        if (Delta.xu[i] < 0 && model.has_ub(i)){
            alpha_primal = std::min(alpha_primal, - It.xu[i] / Delta.xu[i]);
        }
    }
    alpha_primal *= alpha_interior_scaling;

    alpha_dual = 1.0;
    for (int i=0; i<n; ++i){
        if (Delta.zl[i] < 0 && model.has_lb(i)){
            alpha_dual = std::min(alpha_dual, - It.zl[i] / Delta.zl[i]);
        }
        if (Delta.zu[i] < 0 && model.has_ub(i)){
            alpha_dual = std::min(alpha_dual, - It.zu[i] / Delta.zu[i]);
        }
    }
    alpha_dual *= alpha_interior_scaling;

}

