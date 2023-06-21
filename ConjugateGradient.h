#ifndef CONJUGATEGRADIENT_H
#define CONJUGATEGRADIENT_H

#include "SparseMatrix.h"
#include "VectorOperations.h"
#include "NormalEquations.h"

void CG_solve(
    const NormalEquations& A,
    const std::vector<double>& rhs,
    double tol,
    int maxit,
    std::vector<double>& lhs
);




#endif