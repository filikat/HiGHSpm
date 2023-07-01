#ifndef DIRECT_H
#define DIRECT_H

#include "SparseMatrix.h"

void newtonSolve(const HighsSparseMatrix &highs_a,
		 const std::vector<double> &theta,
                 const std::vector<double> &rhs,
		 std::vector<double> &lhs);

/*
void augmentedSolve(const HighsSparseMatrix &highs_a,
                 const std::vector<double> &scaling,
                 const std::vector<double> &rhs,
                 std::vector<double> &lhs);
*/
#endif
