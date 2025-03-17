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

class IpmModel {
 private:
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

  // Put the model into correct formulation
  void reformulate();

  // Scale the problem
  void scale();

 public:
  // Initialize the model
  void init(const int num_var, const int num_con, const double* obj,
            const double* rhs, const double* lower, const double* upper,
            const int* A_ptr, const int* A_rows, const double* A_vals,
            const char* constraints, const std::string& pb_name);

  // Compute range of coefficients
  void checkCoefficients() const;

  // Unscale a given solution
  void unscale(std::vector<double>& x, std::vector<double>& xl,
               std::vector<double>& xu, std::vector<double>& slack,
               std::vector<double>& y, std::vector<double>& zl,
               std::vector<double>& zu) const;
  void unscale(std::vector<double>& x, std::vector<double>& slack,
               std::vector<double>& y, std::vector<double>& z) const;

  double normScaledRhs() const;
  double normScaledObj() const;
  double normUnscaledRhs() const;
  double normUnscaledObj() const;

  // Check if variable has finite lower/upper bound
  bool hasLb(int j) const { return lower_[j] != -kInf; }
  bool hasUb(int j) const { return upper_[j] != kInf; }

  int m() const { return m_; }
  int n() const { return n_; }
  int n_orig() const { return num_var_; }
  const HighsSparseMatrix& A() const { return A_; }
  const std::vector<double>& b() const { return b_; }
  const std::vector<double>& c() const { return c_; }
  double lb(int i) const { return lower_[i]; }
  double ub(int i) const { return upper_[i]; }
  const std::string& name() const { return pb_name_; }
  char constraint(int i) const { return constraints_[i]; }
  int colexp(int i) const { return colexp_[i]; }
  int rowexp(int i) const { return rowexp_[i]; }
  bool ready() const { return ready_; }
  bool scaled() const { return colexp_.size() > 0; }
};

#endif
