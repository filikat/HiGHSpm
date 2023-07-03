#ifndef CONJUGATEGRADIENT_H
#define CONJUGATEGRADIENT_H

#include "NormalEquations.h"
#include "VectorOperations.h"
#include "util/HighsSparseMatrix.h"

void CG_solve(const NormalEquations &A, const std::vector<double> &rhs,
              double tol, int maxit, std::vector<double> &lhs, int *cg_iter);

#endif
