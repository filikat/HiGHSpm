#ifndef DIRECT_H
#define DIRECT_H

#include "experimentData.h"
#include "util/HighsSparseMatrix.h"

void productAThetaAT(const HighsSparseMatrix& matrix, const double* theta,
		     const std::vector<double>& x, std::vector<double>& y);
HighsSparseMatrix computeAThetaAT_inner_product(const HighsSparseMatrix& matrix, const double* theta);
std::ostream& operator<<(std::ostream& os, const ExperimentData& data);

int newtonSolve(const HighsSparseMatrix &highs_a,
		const std::vector<double> &theta,
		const std::vector<double> &rhs,
		std::vector<double> &lhs,
		const int option_max_dense_col,
		const double option_dense_col_tolerance,
		ExperimentData& data
		);

/*
void augmentedSolve(const HighsSparseMatrix &highs_a,
                 const std::vector<double> &scaling,
                 const std::vector<double> &rhs,
                 std::vector<double> &lhs);
*/
#endif
