#include <cassert>
#include <iostream>
#include "Highs.h"
#include "IPM_caller.h"



int main(int argc, char** argv){

    if (argc < 2 || argc > 3){
        std::cerr << "======= How to use: ./ipm LP_name.mps =======\n";
        return 1;
    }


    // ===================================================================================
    // READ PROBLEM
    // ===================================================================================

    // Read LP using Highs MPS read
    Highs highs;
    HighsStatus status = highs.readModel(argv[1]);
    assert(status == HighsStatus::kOk);
    HighsLp lp = highs.getLp();



    // ===================================================================================
    // CHANGE FORMULATION
    // ===================================================================================
    // Input problem must be in the form
    //
    //  min   obj^T * x
    //  s.t.  A * x {<=,=,>=} rhs
    //        lower <= x <= upper
    //
    //  constraints[i] is : 0 for  =
    //                      1 for >=
    //                     -1 for <=
    // ===================================================================================

    // Make a local copy of LP data to be modified
    int n = lp.num_col_;
    int m = lp.num_row_;
    std::vector<double> obj(lp.col_cost_);
    std::vector<double> lower(lp.col_lower_);
    std::vector<double> upper(lp.col_upper_);
    std::vector<int> colptr(lp.a_matrix_.start_);
    std::vector<int> rowind(lp.a_matrix_.index_);
    std::vector<double> values(lp.a_matrix_.value_);

    // Prepare vectors for different formulation
    std::vector<double> rhs(lp.num_row_);
    std::vector<int> constraints(lp.num_row_);
    
    int num_slacks{};

    // Set up constraints and rhs based on row_lower_ and row_upper_ for each constraint
    for (int i=0; i<m; ++i){

        // equality constraint
        if (lp.row_lower_[i] == lp.row_upper_[i]){
            constraints[i] = 0;
            rhs[i] = lp.row_lower_[i];
        }

        // constraint <=
        else if(lp.row_lower_[i] <= -kHighsInf && lp.row_upper_[i] < kHighsIInf){
            constraints[i] = -1;
            rhs[i] = lp.row_upper_[i];
        }

        // constraint >=
        else if(lp.row_lower_[i] > -kHighsInf && lp.row_upper_[i] >= kHighsIInf){
            constraints[i] = 1;
            rhs[i] = lp.row_lower_[i];
        }

        // no constraint
        else if (lp.row_lower_[i] <= -kHighsInf && lp.row_upper_[i] >= kHighsIInf){
            std::cout<<"======= Free variable not yet supported =======\n";
            return 1;
        }

        // general constraint
        else{
            // keep track of how many slacks are needed to allocate memory
            ++num_slacks;
        }

    }

    // Reserve memory for slacks
    obj.reserve(n + num_slacks);
    lower.reserve(n + num_slacks);
    upper.reserve(n + num_slacks);
    colptr.reserve(n + num_slacks + 1);
    rowind.reserve(lp.a_matrix_.numNz() + num_slacks);
    values.reserve(lp.a_matrix_.numNz() + num_slacks);

    for (int i=0; i<m; ++i){

        if ( lp.row_lower_[i] != lp.row_upper_[i] && lp.row_lower_[i] > -kHighsIInf && lp.row_upper_[i] < kHighsIInf){
            // If both row_lower_ and row_upper_ are finite and different, add slack:
            //   lb <= a^T * x <= ub
            //   becomes
            //   a^T * x - s = 0 and lb <= s <= ub.
            //
            //  This requires:
            //   - updating obj, lower, upper
            //   - adding a column of -identity to A

            constraints[i] = 0;
            rhs[i] = 0.0;

            // add slack
            ++n;
            obj.push_back(0.0);
            lower.push_back(lp.row_lower_[i]);
            upper.push_back(lp.row_upper_[i]);
            colptr.push_back(colptr.back()+1);
            rowind.push_back(i);
            values.push_back(-1.0);
        }

    }


    
    // ===================================================================================
    // LOAD AND SOLVE THE PROBLEM
    // ===================================================================================

    // create instance of IPM
    IPM_caller ipm{};

    // load the problem
    ipm.Load(
        n,
        m,
        obj.data(),         
        rhs.data(),
        lower.data(),
        upper.data(),
        colptr.data(),
        rowind.data(),
        values.data(),
        constraints.data()
        );

    // solve LP
    ipm.Solve(100,1e-6);


    


    return 0;
}