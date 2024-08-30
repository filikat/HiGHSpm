#include "IpmModel.h"

// =======================================================================
// RESIZE THE MODEL
// =======================================================================
void IpmModel::resize(int num_var, int num_con) {
  num_var_ = num_var;
  num_con_ = num_con;
  obj_.resize(num_var_, 0.0);
  rhs_.resize(num_con_);
  lower_.resize(num_var_);
  upper_.resize(num_var_);
}
