#ifndef IPM_MODEL_H
#define IPM_MODEL_H

#include <limits>
#include <string>
#include <vector>

#include "CurtisReidScaling.h"
#include "Ipm_aux.h"
#include "util/HighsSparseMatrix.h"

// Optimization problem:
// min  c^T * x
// s.t. A * x = b
//      x - xl = lower
//      x + xu = upper
//      xl, xu >= 0
//
// A is of size num_con x num_var, stored in CSC format using ptr, rows, vals.
//
// See Schork, Gondzio "Implementation of an interior point method with basis
// preconditioning", Math. Prog. Comput. 12, 2020
//

// double infinity
const double kInf = std::numeric_limits<double>::infinity();

struct IpmModel {
  // data of original problem
  int num_var_{};

  // data of reformulated problem
  int n_{};
  int m_{};
  std::vector<double> c_{};
  std::vector<double> b_{};
  std::vector<double> lower_{};
  std::vector<double> upper_{};
  HighsSparseMatrix A_{};
  std::vector<char> constraints_{};
  std::string pb_name_{};

  bool ready_ = false;

  // exponents for scaling
  std::vector<int> colexp_{};
  std::vector<int> rowexp_{};
  int cexp_{};
  int bexp_{};

  // Check if variable has finite lower/upper bound
  bool hasLb(int j) const { return lower_[j] != -kInf; }
  bool hasUb(int j) const { return upper_[j] != kInf; }

  // Initialize the model
  void init(const int num_var, const int num_con, const double* obj,
            const double* rhs, const double* lower, const double* upper,
            const int* A_ptr, const int* A_rows, const double* A_vals,
            const char* constraints, const std::string& pb_name);

  // Put the model into correct formulation
  void reformulate();

  // Compute range of coefficients
  void checkCoefficients() const;

  // (Un)scale the matrix
  void scale();
  void unscale(std::vector<double>& x, std::vector<double>& xl,
               std::vector<double>& xu, std::vector<double>& slack,
               std::vector<double>& y, std::vector<double>& zl,
               std::vector<double>& zu) const;

  double normScaledRhs() const;
  double normScaledObj() const;
  double normUnscaledRhs() const;
  double normUnscaledObj() const;

  friend class Ipm;
};

#endif
