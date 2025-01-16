#include "IpmModel.h"

void IpmModel::init(const int num_var, const int num_con, const double* obj,
                    const double* rhs, const double* lower, const double* upper,
                    const int* A_ptr, const int* A_rows, const double* A_vals,
                    const int* constraints, const std::string& pb_name) {
  // copy the input into the model

  num_var_ = num_var;
  num_con_ = num_con;
  c_ = std::vector<double>(obj, obj + num_var_);
  b_ = std::vector<double>(rhs, rhs + num_con_);
  lower_ = std::vector<double>(lower, lower + num_var_);
  upper_ = std::vector<double>(upper, upper + num_var_);

  int Annz = A_ptr[num_var_];
  A_.num_col_ = num_var_;
  A_.num_row_ = num_con_;
  A_.start_ = std::vector<int>(A_ptr, A_ptr + num_var_ + 1);
  A_.index_ = std::vector<int>(A_rows, A_rows + Annz);
  A_.value_ = std::vector<double>(A_vals, A_vals + Annz);

  constraints_ = std::vector<int>(constraints, constraints + num_con_);

  pb_name_ = pb_name;

  ready_ = true;
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
      c_.push_back(0.0);

      // add column of identity to A_
      std::vector<int> temp_ind{i};
      std::vector<double> temp_val{1.0};
      A_.addVec(1, temp_ind.data(), temp_val.data());

      // set scaling to 1, i.e. exponent to zero
      if (colexp_.size() > 0) {
        colexp_.push_back(0);
      }
    }
  }
}

void IpmModel::checkCoefficients() const {
  // compute max and min entry of A in absolute value
  double Amin = kInf;
  double Amax = 0.0;
  for (int col = 0; col < A_.num_col_; ++col) {
    for (int el = A_.start_[col]; el < A_.start_[col + 1]; ++el) {
      double val = std::abs(A_.value_[el]);
      if (val != 0.0) {
        Amin = std::min(Amin, val);
        Amax = std::max(Amax, val);
      }
    }
  }
  if (Amin == kInf) Amin = 0.0;

  // compute max and min entry of c
  double cmin = kInf;
  double cmax = 0.0;
  for (int i = 0; i < num_var_; ++i) {
    if (c_[i] != 0.0) {
      cmin = std::min(cmin, std::abs(c_[i]));
      cmax = std::max(cmax, std::abs(c_[i]));
    }
  }
  if (cmin == kInf) cmin = 0.0;

  // compute max and min entry of b
  double bmin = kInf;
  double bmax = 0.0;
  for (int i = 0; i < num_con_; ++i) {
    if (b_[i] != 0.0) {
      bmin = std::min(bmin, std::abs(b_[i]));
      bmax = std::max(bmax, std::abs(b_[i]));
    }
  }
  if (bmin == kInf) bmin = 0.0;

  // compute max and min for bounds
  double boundmin = kInf;
  double boundmax = 0.0;
  for (int i = 0; i < num_var_; ++i) {
    /*if (std::isfinite(lower_[i]) && std::isfinite(upper_[i]) &&
        upper_[i] != lower_[i]) {
      boundmin = std::min(boundmin, std::abs(upper_[i] - lower_[i]));
      boundmax = std::max(boundmax, std::abs(upper_[i] - lower_[i]));
    }*/

    if (lower_[i] != 0.0 && std::isfinite(lower_[i])) {
      boundmin = std::min(boundmin, std::abs(lower_[i]));
      boundmax = std::max(boundmax, std::abs(lower_[i]));
    }
    if (upper_[i] != 0.0 && std::isfinite(upper_[i])) {
      boundmin = std::min(boundmin, std::abs(upper_[i]));
      boundmax = std::max(boundmax, std::abs(upper_[i]));
    }
  }
  if (boundmin == kInf) boundmin = 0.0;

  // compute max and min scaling
  double scalemin = kInf;
  double scalemax = 0.0;
  if (colexp_.size() > 0) {
    for (int i = 0; i < num_var_; ++i) {
      double scaling = std::ldexp(1.0, colexp_[i]);
      scalemin = std::min(scalemin, scaling);
      scalemax = std::max(scalemax, scaling);
    }
  }
  if (rowexp_.size() > 0) {
    for (int i = 0; i < num_con_; ++i) {
      double scaling = std::ldexp(1.0, rowexp_[i]);
      scalemin = std::min(scalemin, scaling);
      scalemax = std::max(scalemax, scaling);
    }
  }
  if (scalemin == kInf) scalemin = 0.0;

  // print ranges
  printf("Range of A      : [%5.1e, %5.1e], ratio ", Amin, Amax);
  Amin == 0.0 ? printf("-\n") : printf("%.1e\n", Amax / Amin);
  printf("Range of b      : [%5.1e, %5.1e], ratio ", bmin, bmax);
  bmin == 0.0 ? printf("-\n") : printf("%.1e\n", bmax / bmin);
  printf("Range of c      : [%5.1e, %5.1e], ratio ", cmin, cmax);
  cmin == 0.0 ? printf("-\n") : printf("%.1e\n", cmax / cmin);
  printf("Range of bounds : [%5.1e, %5.1e], ratio ", boundmin, boundmax);
  boundmin == 0.0 ? printf("-\n") : printf("%.1e\n", boundmax / boundmin);
  printf("Scaling coeff   : [%5.1e, %5.1e], ratio ", scalemin, scalemax);
  scalemin == 0.0 ? printf("-\n") : printf("%.1e\n", scalemax / scalemin);
}

void IpmModel::scale() {
  // Apply Curtis-Reid scaling and scale the problem accordingly

  // *********************************************************************
  // Compute scaling
  // *********************************************************************
  // Transformation:
  // A -> R * A * C
  // b -> beta * R * b
  // c -> gamma * C * c
  // x -> beta * C^-1 * x
  // y -> gamma * R^-1 * y
  // z -> gamma * C * z
  // where R is row scaling, C is col scaling, beta is unif scaling of b, gamma
  // is unif scaling of c.

  // compute exponents for CR scaling of matrix A
  colexp_.resize(num_var_);
  rowexp_.resize(num_con_);
  CurtisReidScaling(A_.start_, A_.index_, A_.value_, rowexp_, colexp_);

  // compute uniform scaling of b and c
  // (it doesn't seem to help and it can actually hamper convergence)
  /*
  double bmax = 0.0;
  double bmin = std::numeric_limits<double>::max();
  for (int i = 0; i < num_con_; ++i) {
    if (b_[i] != 0.0) {
      bmax = std::max(bmax, std::abs(b_[i]));
      bmin = std::min(bmin, std::abs(b_[i]));
    }
  }
  std::frexp(1.0 / sqrt(bmax * bmin), &bexp_);
  double cmax = 0.0;
  double cmin = std::numeric_limits<double>::max();
  for (int i = 0; i < num_var_; ++i) {
    if (c_[i] != 0.0) {
      cmax = std::max(cmax, std::abs(c_[i]));
      cmin = std::min(cmin, std::abs(c_[i]));
    }
  }
  std::frexp(1.0 / sqrt(cmax * cmin), &cexp_);
  */

  // *********************************************************************
  // Apply scaling
  // *********************************************************************
  // The scaling is given by exponents.
  // To multiply by the scaling a quantity x: std::ldexp(x,  exp).
  // To divide   by the scaling a quantity x: std::ldexp(x, -exp).
  // Using ldexp instead of * or / ensures that only the exponent bits of the
  // floating point number are manipulated.
  //
  // Each row and each columns are scaled by their own exponent.
  // There are uniform scalings of c_ and b_.

  // Column has been scaled up by colscale_[col], so cost is scaled up and
  // bounds are scaled down
  for (int col = 0; col < num_var_; ++col) {
    c_[col] = std::ldexp(c_[col], colexp_[col]);
    c_[col] = std::ldexp(c_[col], cexp_);

    lower_[col] = std::ldexp(lower_[col], -colexp_[col]);
    lower_[col] = std::ldexp(lower_[col], bexp_);

    upper_[col] = std::ldexp(upper_[col], -colexp_[col]);
    upper_[col] = std::ldexp(upper_[col], bexp_);
  }

  // Row has been scaled up by rowscale_[row], so b is scaled up
  for (int row = 0; row < num_con_; ++row) {
    b_[row] = std::ldexp(b_[row], rowexp_[row]);
    b_[row] = std::ldexp(b_[row], bexp_);
  }

  // Each entry of the matrix is scaled by the corresponding row and col factor
  for (int col = 0; col < num_var_; ++col) {
    for (int el = A_.start_[col]; el < A_.start_[col + 1]; ++el) {
      int row = A_.index_[el];
      A_.value_[el] = std::ldexp(A_.value_[el], rowexp_[row]);
      A_.value_[el] = std::ldexp(A_.value_[el], colexp_[col]);
    }
  }
}

void IpmModel::unscale(Iterate& it) {
  // Undo the scaling

  if (colexp_.size() > 0) {
    for (int i = 0; i < num_var_; ++i) {
      it.x[i] = std::ldexp(it.x[i], colexp_[i]);
      it.x[i] = std::ldexp(it.x[i], -bexp_);

      it.xl[i] = std::ldexp(it.xl[i], colexp_[i]);
      it.xl[i] = std::ldexp(it.xl[i], -bexp_);

      it.xu[i] = std::ldexp(it.xu[i], colexp_[i]);
      it.xu[i] = std::ldexp(it.xu[i], -bexp_);

      it.zl[i] = std::ldexp(it.zl[i], -colexp_[i]);
      it.zl[i] = std::ldexp(it.zl[i], -cexp_);

      it.zu[i] = std::ldexp(it.zu[i], -colexp_[i]);
      it.zu[i] = std::ldexp(it.zu[i], -cexp_);
    }
  }
  if (rowexp_.size() > 0) {
    for (int i = 0; i < num_con_; ++i) {
      it.y[i] = std::ldexp(it.y[i], rowexp_[i]);
      it.y[i] = std::ldexp(it.y[i], -cexp_);
    }
  }

  // set variables that were ignored
  for (int i = 0; i < num_var_; ++i) {
    if (!hasLb(i)) {
      it.xl[i] = kInf;
      it.zl[i] = kInf;
    }
    if (!hasUb(i)) {
      it.xu[i] = kInf;
      it.zu[i] = kInf;
    }
  }
}

double IpmModel::normScaledRhs() const {
  double norm_rhs = infNorm(b_);
  for (double d : lower_)
    if (std::isfinite(d)) norm_rhs = std::max(norm_rhs, std::abs(d));
  for (double d : upper_)
    if (std::isfinite(d)) norm_rhs = std::max(norm_rhs, std::abs(d));
  return norm_rhs;
}

double IpmModel::normScaledObj() const { return infNorm(c_); }

double IpmModel::normUnscaledObj() const {
  double norm_obj = 0.0;
  for (int i = 0; i < num_var_; ++i) {
    double val = std::abs(c_[i]);
    if (colexp_.size() > 0) val = std::ldexp(val, -cexp_ - colexp_[i]);
    norm_obj = std::max(norm_obj, val);
  }
  return norm_obj;
}

double IpmModel::normUnscaledRhs() const {
  double norm_rhs = 0.0;
  for (int i = 0; i < num_con_; ++i) {
    double val = std::abs(b_[i]);
    if (rowexp_.size() > 0) val = std::ldexp(val, -bexp_ - rowexp_[i]);
    norm_rhs = std::max(norm_rhs, val);
  }
  for (int i = 0; i < num_var_; ++i) {
    if (std::isfinite(lower_[i])) {
      double val = std::abs(lower_[i]);
      if (colexp_.size() > 0) val = std::ldexp(val, colexp_[i] - bexp_);
      norm_rhs = std::max(norm_rhs, val);
    }
    if (std::isfinite(upper_[i])) {
      double val = std::abs(upper_[i]);
      if (colexp_.size() > 0) val = std::ldexp(val, colexp_[i] - bexp_);
      norm_rhs = std::max(norm_rhs, val);
    }
  }
  return norm_rhs;
}
