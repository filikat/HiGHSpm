#ifndef IPM_MODEL_H
#define IPM_MODEL_H

#include <limits>
#include <string>
#include <vector>

#include "util/HighsSparseMatrix.h"

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
const double kInf = std::numeric_limits<double>::infinity();

// Constraint types
const double kConstraintTypeLower = 1;
const double kConstraintTypeEqual = 0;
const double kConstraintTypeUpper = -1;

class IpmModel {
  int num_var_{};
  int num_con_{};
  std::vector<double> obj_{};
  std::vector<double> rhs_{};
  std::vector<double> lower_{};
  std::vector<double> upper_{};
  HighsSparseMatrix A_{};
  std::string pb_name_{};

 public:
  // =======================================================================
  // RESIZE THE MODEL
  // =======================================================================
  // Allocate space for model with given number of variables and constraints
  void resize(int num_var,  // number of variables in the model
              int num_con   // number of constraints in the model
  );

  // =======================================================================
  // FIND FINITE BOUNDS
  // =======================================================================
  // check if variable has finite lower/upper bound
  bool hasLb(int j) const { return lower_[j] != -kInf; }
  bool hasUb(int j) const { return upper_[j] != kInf; }

  friend class Ipm;
};

#endif
