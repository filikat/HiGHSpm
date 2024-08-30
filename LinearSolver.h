#ifndef LINEAR_SOLVER_H
#define LINEAR_SOLVER_H

#include "Ipm_aux.h"
#include "Ipm_const.h"
#include "VectorOperations.h"
#include "util/HighsSparseMatrix.h"
#include <vector>

// Interface class for solving augmented system or normal equations.
// Any linear solver needs to define the functions:
// - factorAS: to factorize the augmented system
// - solveAS: to solve a linear system with the augmented system
// - factorNE: to factorize the normal equations
// - solveNE: to solve a linear system with the normal equations
// - clear: to reset the data structure for the next factorization.

class LinearSolver {
public:
  bool valid_ = false;

  virtual int factorAS(const HighsSparseMatrix &A,
                       const std::vector<double> &theta) = 0;

  virtual int
  solveAS(const HighsSparseMatrix &A, const std::vector<double> &theta,
          const std::vector<double> &rhs_x, const std::vector<double> &rhs_y,
          std::vector<double> &lhs_x, std::vector<double> &lhs_y) = 0;

  virtual int factorNE(const HighsSparseMatrix &A,
                       const std::vector<double> &theta) = 0;

  virtual int solveNE(const HighsSparseMatrix &A,
                      const std::vector<double> &theta,
                      const std::vector<double> &rhs,
                      std::vector<double> &lhs) = 0;

  virtual void clear() = 0;
};

#endif