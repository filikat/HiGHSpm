#include "IPM_model.h"



// =======================================================================
// RESIZE THE MODEL
// =======================================================================
void IPM_model::resize(int input_var,int input_con){
    num_var = input_var;
    num_con = input_con;
    obj.resize(num_var,0.0);
    rhs.resize(num_con);
    lower.resize(num_var);
    upper.resize(num_var);
}


