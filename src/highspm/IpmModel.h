#ifndef HIGHSPM_IPM_MODEL_H
#define HIGHSPM_IPM_MODEL_H

#include <limits>
#include <string>
#include <vector>

#include "CurtisReidScaling.h"
#include "ipm/ipx/lp_solver.h"
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

class IpmModel {
 private:
  // data of original problem
  int n_orig_{};
  int m_orig_{};
  const double* c_orig_;
  const double* b_orig_;
  const double* lower_orig_;
  const double* upper_orig_;
  const int* A_ptr_orig_;
  const int* A_rows_orig_;
  const double* A_vals_orig_;
  const char* constraints_orig_;
  double offset_;

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

  // coefficients for scaling
  std::vector<double> colscale_{};
  std::vector<double> rowscale_{};

  // information about empty rows, for postprocessing
  std::vector<int> rows_shift_{};

  // Put the model into correct formulation
  void reformulate();

  // Scale the problem
  void scale();

  void preprocess();

 public:
  // Initialize the model
  void init(const int num_var, const int num_con, const double* obj,
            const double* rhs, const double* lower, const double* upper,
            const int* A_ptr, const int* A_rows, const double* A_vals,
            const char* constraints, double offset, const std::string& pb_name);

  // Compute range of coefficients
  void checkCoefficients() const;

  void postprocess(std::vector<double>& slack, std::vector<double>& y) const;

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
  bool hasLb(int j) const { return std::isfinite(lower_[j]); }
  bool hasUb(int j) const { return std::isfinite(upper_[j]); }

  int m() const { return m_; }
  int n() const { return n_; }
  int n_orig() const { return n_orig_; }
  const HighsSparseMatrix& A() const { return A_; }
  const std::vector<double>& b() const { return b_; }
  const std::vector<double>& c() const { return c_; }
  double lb(int i) const { return lower_[i]; }
  double ub(int i) const { return upper_[i]; }
  const std::string& name() const { return pb_name_; }
  char constraint(int i) const { return constraints_[i]; }
  double colScale(int i) const { return colscale_[i]; }
  double rowScale(int i) const { return rowscale_[i]; }
  bool ready() const { return ready_; }
  bool scaled() const { return colscale_.size() > 0; }
  double offset() const { return offset_; }

  int loadIntoIpx(ipx::LpSolver& lps) const;
};

#endif
