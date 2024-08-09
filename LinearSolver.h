#ifndef LINEAR_SOLVER_H
#define LINEAR_SOLVER_H

#include "IPM_aux.h"
#include "IPM_const.h"
#include "VectorOperations.h"
#include "util/HighsSparseMatrix.h"
#include <vector>

// Interface class for solving augmented system or normal equations.
// Any linear solver needs to define the functions:
// - FactorAS: to factorize the augmented system
// - SolveAS: to solve a linear system with the augmented system
// - FactorNE: to factorize the normal equations
// - SolveNE: to solve a linear system with the normal equations
// - Clear: to reset the data structure for the next factorization.

class LinearSolver {
public:
  bool valid = false;

  virtual int FactorAS(const HighsSparseMatrix &A,
                       const std::vector<double> &theta) = 0;

  virtual int
  SolveAS(const HighsSparseMatrix &A, const std::vector<double> &theta,
          const std::vector<double> &rhs_x, const std::vector<double> &rhs_y,
          std::vector<double> &lhs_x, std::vector<double> &lhs_y) = 0;

  virtual int FactorNE(const HighsSparseMatrix &A,
                       const std::vector<double> &theta) = 0;

  virtual int SolveNE(const HighsSparseMatrix &A,
                      const std::vector<double> &theta,
                      const std::vector<double> &rhs,
                      std::vector<double> &lhs) = 0;

  virtual void Clear() = 0;
};

#endif