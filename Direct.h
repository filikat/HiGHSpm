#ifndef DIRECT_H
#define DIRECT_H

#include "ExperimentData.h"
#include "util/HighsSparseMatrix.h"

bool increasing_index(const HighsSparseMatrix& matrix);
void productAThetaAT(const HighsSparseMatrix& matrix,
		     const std::vector<double>& theta,
		     const std::vector<double>& x,
		     std::vector<double>& result);
HighsSparseMatrix computeAThetaAT(const HighsSparseMatrix& matrix,
				  const std::vector<double>& theta);

int augmentedSolve(const HighsSparseMatrix &highs_a,
		    const std::vector<double> &scaling,
		    const std::vector<double> &rhs_x,
		    const std::vector<double> &rhs_y,
		    std::vector<double> &lhs_x,
		    std::vector<double> &lhs_y,
		    ExperimentData& data);

int newtonSolve(const HighsSparseMatrix &highs_a,
		const std::vector<double> &theta,
		const std::vector<double> &rhs,
		std::vector<double> &lhs,
		const int option_max_dense_col,
		const double option_dense_col_tolerance,
		ExperimentData& data);

#endif
