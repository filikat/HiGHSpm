#ifndef IPM_MODEL_H
#define IPM_MODEL_H

#include "util/HighsSparseMatrix.h"
#include <limits>
#include <vector>
#include <string>

// Optimization problem:
// min obj^T * x
// s.t. A * x = rhs
//      x - xl = lower
//      x + xu = upper
//      xl, xu >= 0
//
// A is of size num_con x num_var and
// is stored in CSC format using:
// colptr, rowind, values
//
// See Schork, Gondzio "Implementation of an interior point method with basis
// preconditioning", Math. Prog. Comput. 12, 2020
//

// double infinity
const double INF = std::numeric_limits<double>::infinity();

// Constraint types
const double kConstraintTypeLower = 1;
const double kConstraintTypeEqual = 0;
const double kConstraintTypeUpper = -1;

class IPM_model {

  int num_var{};
  int num_con{};
  std::vector<double> obj{};
  std::vector<double> rhs{};
  std::vector<double> lower{};
  std::vector<double> upper{};
  HighsSparseMatrix highs_a{};
  std::string pb_name{};

public:
  // =======================================================================
  // RESIZE THE MODEL
  // =======================================================================
  // Allocate space for model with given number of variables and constraints
  void resize(int input_var, // number of variables in the model
              int input_con  // number of constraints in the model
  );

  // =======================================================================
  // FIND FINITE BOUNDS
  // =======================================================================
  // check if variable has finite lower/upper bound
  bool has_lb(int j) const { return lower[j] != -INF; }
  bool has_ub(int j) const { return upper[j] != INF; }

  friend class IPM_caller;
};

#endif
