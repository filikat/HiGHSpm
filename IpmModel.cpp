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

      // set scaling to 1
      if (colscale_.size() > 0) {
        colscale_.push_back(1.0);
      }
    }
  }
}

void IpmModel::checkCoefficients() {
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

  // compute max and min entry of obj
  double cmin = kInf;
  double cmax = 0.0;
  for (int i = 0; i < num_var_; ++i) {
    if (obj_[i] != 0.0) {
      cmin = std::min(cmin, std::abs(obj_[i]));
      cmax = std::max(cmax, std::abs(obj_[i]));
    }
  }
  if (cmin == kInf) cmin = 0.0;

  // compute max and min entry of rhs
  double bmin = kInf;
  double bmax = 0.0;
  for (int i = 0; i < num_con_; ++i) {
    if (rhs_[i] != 0.0) {
      bmin = std::min(bmin, std::abs(rhs_[i]));
      bmax = std::max(bmax, std::abs(rhs_[i]));
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
  for (int i = 0; i < num_var_; ++i) {
    scalemin = std::min(scalemin, colscale_[i]);
    scalemax = std::max(scalemax, colscale_[i]);
  }
  for (int i = 0; i < num_con_; ++i) {
    scalemin = std::min(scalemin, rowscale_[i]);
    scalemax = std::max(scalemax, rowscale_[i]);
  }

  // print ranges
  printf("\nCoefficients range\n");
  printf("Range of A     : [%5.1e, %5.1e], ratio %.1e\n", Amin, Amax,
         Amax / Amin);
  printf("Range of b     : [%5.1e, %5.1e], ratio %.1e\n", bmin, bmax,
         bmax / bmin);
  printf("Range of c     : [%5.1e, %5.1e], ratio %.1e\n", cmin, cmax,
         cmax / cmin);
  printf("Range of bounds: [%5.1e, %5.1e], ratio %.1e\n", boundmin, boundmax,
         boundmax / boundmin);
  printf("Scaling coeff  : [%5.1e, %5.1e]\n", scalemin, scalemax);
  printf("\n");
}

void IpmModel::scale() {
  // Apply matrix equilibration and scale the problem accordingly

  if (equilibrate()) {
    // Matrix has been scaled

    for (int col = 0; col < num_var_; ++col) {
      // Column has been scaled up by colscale_[col], so cost is scaled up and
      // bounds are scaled down
      obj_[col] *= colscale_[col];
      lower_[col] /= colscale_[col];
      upper_[col] /= colscale_[col];
    }

    for (int row = 0; row < num_con_; ++row) {
      // Row has been scaled up by rowscale_[row], so rhs is scaled up
      rhs_[row] *= rowscale_[row];
    }
  }
}

void IpmModel::scaleCR() {
  // Apply Curtis-Reid scaling and scale the problem accordingly

  double objscale{};
  double rhsscale{};
  colscale_.resize(num_var_);
  rowscale_.resize(num_con_);

  CurtisReidScaling(A_.start_, A_.index_, A_.value_, rhs_, obj_, objscale,
                    rhsscale, rowscale_, colscale_);

  // Column has been scaled up by colscale_[col], so cost is scaled up and
  // bounds are scaled down
  for (int col = 0; col < num_var_; ++col) {
    obj_[col] *= colscale_[col];
    lower_[col] /= colscale_[col];
    upper_[col] /= colscale_[col];
  }

  // Row has been scaled up by rowscale_[row], so rhs is scaled up
  for (int row = 0; row < num_con_; ++row) {
    rhs_[row] *= rowscale_[row];
  }

  // Each entry of the matrix is scaled by the corresponding row and col factor
  for (int col = 0; col < num_var_; ++col) {
    for (int el = A_.start_[col]; el < A_.start_[col + 1]; ++el) {
      int row = A_.index_[el];
      A_.value_[el] *= rowscale_[row];
      A_.value_[el] *= colscale_[col];
    }
  }
}

void IpmModel::unscale(Iterate& it) {
  // Undo the scaling

  if (colscale_.size() > 0) {
    for (int i = 0; i < num_var_; ++i) {
      it.x[i] *= colscale_[i];
      it.xl[i] *= colscale_[i];
      it.xu[i] *= colscale_[i];
      it.zl[i] /= colscale_[i];
      it.zu[i] /= colscale_[i];
    }
  }
  if (rowscale_.size() > 0) {
    for (int i = 0; i < num_con_; ++i) {
      it.y[i] *= rowscale_[i];
    }
  }
}

double equilFactor(int exp, int min, int max) {
  // Equilibration factor needed by matrix equilibration
  if (exp < min) {
    return std::ldexp(1.0, (min - exp + 1) / 2);
  }
  if (exp > max) {
    return std::ldexp(1.0, -((exp - max + 1) / 2));
  }
  return 1.0;
}

bool IpmModel::equilibrate() {
  // Iterative procedure to get better scaled matrix A:
  // ideally, each entry of A should be written as x * 2^exp,
  // with x \in [0.5,1.0) and exp \in [expmin,expmax].
  // Taken from IPX Model::EquilibrateMatrix().

  const int expmin = 0;
  const int expmax = 3;
  const int maxiter = 10;

  bool out_of_range = false;

  colscale_.resize(0);
  rowscale_.resize(0);

  int exp;

  // check if matrix needs equilibration
  for (int el = 0; el < A_.numNz(); ++el) {
    std::frexp(std::abs(A_.value_[el]), &exp);
    if (exp < expmin || exp > expmax) {
      out_of_range = true;
      break;
    }
  }
  if (!out_of_range) return false;

  // initialize equilibration factors
  colscale_.resize(num_var_, 1.0);
  rowscale_.resize(num_con_, 1.0);

  // initialize vectors for max entry in rows and cols
  std::vector<double> colmax(num_var_);
  std::vector<double> rowmax(num_con_);

  // iterate
  for (int iter = 0; iter < maxiter; ++iter) {
    // compute \infty norm of rows and cols
    rowmax.assign(num_con_, 0.0);
    colmax.assign(num_var_, 0.0);
    for (int col = 0; col < num_var_; ++col) {
      for (int el = A_.start_[col]; el < A_.start_[col + 1]; ++el) {
        int row = A_.index_[el];
        double val = std::abs(A_.value_[el]);
        colmax[col] = std::max(colmax[col], val);
        rowmax[row] = std::max(rowmax[row], val);
      }
    }

    out_of_range = false;

    int count_rows{};
    int count_cols{};

    // compute scaling factors for rows
    for (int row = 0; row < num_con_; ++row) {
      std::frexp(rowmax[row], &exp);
      rowmax[row] = equilFactor(exp, expmin, expmax);
      if (rowmax[row] != 1.0) {
        out_of_range = true;
        rowscale_[row] *= rowmax[row];
        ++count_rows;
      }
    }

    // compute scaling factors for columns
    for (int col = 0; col < num_var_; ++col) {
      std::frexp(colmax[col], &exp);
      colmax[col] = equilFactor(exp, expmin, expmax);
      if (colmax[col] != 1.0) {
        out_of_range = true;
        colscale_[col] *= colmax[col];
        ++count_cols;
      }
    }

    // matrix has been equilibrated, no need to keep iterating
    if (!out_of_range) break;

    // apply scaling to matrix
    for (int col = 0; col < num_var_; ++col) {
      for (int el = A_.start_[col]; el < A_.start_[col + 1]; ++el) {
        int row = A_.index_[el];
        A_.value_[el] *= colmax[col];
        A_.value_[el] *= rowmax[row];
      }
    }

    printf("Scaled %d rows, %d cols\n", count_rows, count_cols);
  }

  return true;
}