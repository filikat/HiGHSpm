#include "IpmModel.h"

void IpmModel::init(const int num_var, const int num_con, const double* obj,
                    const double* rhs, const double* lower, const double* upper,
                    const int* A_colptr, const int* A_rowind,
                    const double* A_values, const int* constraints,
                    const std::string& pb_name) {
  // copy the input into the model

  num_var_ = num_var;
  num_con_ = num_con;
  obj_ = std::vector<double>(obj, obj + num_var_);
  rhs_ = std::vector<double>(rhs, rhs + num_con_);
  lower_ = std::vector<double>(lower, lower + num_var_);
  upper_ = std::vector<double>(upper, upper + num_var_);

  int Annz = A_colptr[num_var_];
  A_.num_col_ = num_var_;
  A_.num_row_ = num_con_;
  A_.start_ = std::vector<int>(A_colptr, A_colptr + num_var_ + 1);
  A_.index_ = std::vector<int>(A_rowind, A_rowind + Annz);
  A_.value_ = std::vector<double>(A_values, A_values + Annz);

  constraints_ = std::vector<int>(constraints, constraints + num_con_);

  pb_name_ = pb_name;
}

void IpmModel::reformulate() {
  // put the model into correct formulation

  int Annz = A_.numNz();

  for (int i = 0; i < num_con_; ++i) {
    if (constraints_[i] != kConstraintTypeEqual) {
      // inequality constraint, add slack variable

      ++num_var_;

      // lower/upper bound for new slack
      if (constraints_[i] == kConstraintTypeLower) {
        lower_.push_back(-kInf);
        upper_.push_back(0.0);
      } else {
        lower_.push_back(0.0);
        upper_.push_back(kInf);
      }

      // cost for new slack
      obj_.push_back(0.0);

      // add column of identity to A_
      std::vector<int> temp_ind{i};
      std::vector<double> temp_val{1.0};
      A_.addVec(1, temp_ind.data(), temp_val.data());
    }
  }
}

void IpmModel::checkCoefficients() {
  // compute max and min entry of A in absolute value
  double Amin = kInf;
  double Amax = 0.0;
  for (int col = 0; col < A_.num_col_; ++col) {
    for (int el = A_.start_[col]; el < A_.start_[col + 1]; ++el) {
      double val = std::fabs(A_.value_[el]);
      if (val != 0.0) {
        Amin = std::min(Amin, val);
        Amax = std::max(Amax, val);
      }
    }
  }
  if (Amin == kInf) Amin = 0.0;

  // compute max and min entry of obj
  double cmin = kInf;
  double cmax = 0.0;
  for (int i = 0; i < num_var_; ++i) {
    if (obj_[i] != 0.0) {
      cmin = std::min(cmin, std::fabs(obj_[i]));
      cmax = std::max(cmax, std::fabs(obj_[i]));
    }
  }
  if (cmin == kInf) cmin = 0.0;

  // compute max and min entry of rhs
  double bmin = kInf;
  double bmax = 0.0;
  for (int i = 0; i < num_con_; ++i) {
    if (rhs_[i] != 0.0) {
      bmin = std::min(bmin, std::fabs(rhs_[i]));
      bmax = std::max(bmax, std::fabs(rhs_[i]));
    }
  }
  if (bmin == kInf) bmin = 0.0;

  // compute max and min entry of lb,ub
  double boundmin = kInf;
  double boundmax = 0.0;
  for (int i = 0; i < num_var_; ++i) {
    if (lower_[i] != 0.0 && std::isfinite(lower_[i])) {
      boundmin = std::min(boundmin, std::fabs(lower_[i]));
      boundmax = std::max(boundmax, std::fabs(lower_[i]));
    }
    if (upper_[i] != 0.0 && std::isfinite(upper_[i])) {
      boundmin = std::min(boundmin, std::fabs(upper_[i]));
      boundmax = std::max(boundmax, std::fabs(upper_[i]));
    }
  }
  if (boundmin == kInf) boundmin = 0.0;

  // print ranges
  printf("\nCoefficients range\n");
  printf("Range of A     : [%5.1e, %5.1e]\n", Amin, Amax);
  printf("Range of b     : [%5.1e, %5.1e]\n", bmin, bmax);
  printf("Range of c     : [%5.1e, %5.1e]\n", cmin, cmax);
  printf("Range of bounds: [%5.1e, %5.1e]\n", boundmin, boundmax);
  printf("\n");
}