#ifndef HIGHSPM_HPM_MODEL_H
#define HIGHSPM_HPM_MODEL_H

#include <limits>
#include <string>
#include <vector>

#include "CurtisReidScaling.h"
#include "auxiliary/IntConfig.h"
#include "ipm/ipx/lp_solver.h"
#include "util/HighsSparseMatrix.h"

namespace highspm {

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

class HpmModel {
 private:
  // data of original problem
  Int n_orig_{};
  Int m_orig_{};
  const double* c_orig_;
  const double* b_orig_;
  const double* lower_orig_;
  const double* upper_orig_;
  const Int* A_ptr_orig_;
  const Int* A_rows_orig_;
  const double* A_vals_orig_;
  const char* constraints_orig_;
  double offset_;

  // data of reformulated problem
  Int n_{};
  Int m_{};
  std::vector<double> c_{};
  std::vector<double> b_{};
  std::vector<double> lower_{};
  std::vector<double> upper_{};
  HighsSparseMatrix A_{};
  std::vector<char> constraints_{};
  Int num_dense_cols_{};
  double max_col_density_{};

  bool ready_ = false;

  // coefficients for scaling
  std::vector<double> colscale_{};
  std::vector<double> rowscale_{};

  // information about empty rows, for postprocessing
  std::vector<Int> rows_shift_{};
  Int empty_rows_{};

  // Put the model into correct formulation
  void reformulate();

  // Scale the problem
  void scale();

  void preprocess();

  void denseColumns();

 public:
  // Initialise the model
  void init(const Int num_var, const Int num_con, const double* obj,
            const double* rhs, const double* lower, const double* upper,
            const Int* A_ptr, const Int* A_rows, const double* A_vals,
            const char* constraints, double offset);

  // Print information of model
  void print() const;

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
  bool hasLb(Int j) const { return std::isfinite(lower_[j]); }
  bool hasUb(Int j) const { return std::isfinite(upper_[j]); }

  Int m() const { return m_; }
  Int n() const { return n_; }
  Int n_orig() const { return n_orig_; }
  const HighsSparseMatrix& A() const { return A_; }
  const std::vector<double>& b() const { return b_; }
  const std::vector<double>& c() const { return c_; }
  double lb(Int i) const { return lower_[i]; }
  double ub(Int i) const { return upper_[i]; }
  char constraint(Int i) const { return constraints_[i]; }
  double colScale(Int i) const { return colscale_[i]; }
  double rowScale(Int i) const { return rowscale_[i]; }
  bool ready() const { return ready_; }
  bool scaled() const { return colscale_.size() > 0; }
  double offset() const { return offset_; }
  double maxColDensity() const { return max_col_density_; }
  Int numDenseCols() const { return num_dense_cols_; }

  Int loadIntoIpx(ipx::LpSolver& lps) const;
};

}  // namespace highspm

#endif
