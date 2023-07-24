#include "Direct.h"
#include <cmath>

bool increasing_index(const HighsSparseMatrix& matrix) {
  if (matrix.isRowwise()) {
    for (int iRow = 0; iRow < matrix.num_row_; iRow++)
      for (int iEl = matrix.start_[iRow]+1; iEl < matrix.start_[iRow+1]; iRow++) 
	if (matrix.index_[iEl] <= matrix.index_[iEl-1]) return false;
  } else {
    for (int iCol = 0; iCol < matrix.num_col_; iCol++)
      for (int iEl = matrix.start_[iCol]+1; iEl < matrix.start_[iCol+1]; iCol++)
	if (matrix.index_[iEl] <= matrix.index_[iEl-1]) return false;
  }
  return true;
}

HighsSparseMatrix computeAThetaAT_inner_product(const HighsSparseMatrix& matrix, const double* theta){
    HighsSparseMatrix AAT;
    HighsSparseMatrix AT = matrix;
    AT.ensureRowwise();
    assert(increasing_index(AT));
    int AAT_dim = matrix.num_row_;
    AAT.num_col_ = AAT_dim;
    AAT.num_row_ = AAT_dim;
    AAT.start_.resize(AAT_dim+1,0);
 
    std::vector<std::tuple<int, int, double>> non_zero_values;

    // First pass to calculate the number of non-zero elements in each column
    //
    // Value to add for implicit identity matrix
    const double identity_term = 1;
    for (int i = 0; i < AAT_dim; ++i) {
      // If theta is a not a null pointer set theta_i = theta[i], otherwise it's set to 1,
      const double theta_i = theta ? theta[i] : 1;
      for (int j = i; j < AAT_dim; ++j) {
	double dot = 0.0;
	int k = AT.start_[i];
	int l = AT.start_[j];
	while (k < AT.start_[i + 1] && l < AT.start_[j + 1]) {
	  if (AT.index_[k] < AT.index_[l]) {
	    ++k;
	  } else if (AT.index_[k] > AT.index_[l]) {
	    ++l;
	  } else {                 
	    dot += theta_i * AT.value_[k] * AT.value_[l];
	    ++k;
	    ++l;
	  }
	}
	if (dot != 0.0) {
	  non_zero_values.emplace_back(i, j, dot);
	  AAT.start_[i+1]++;
	  if (i != j) AAT.start_[j+1]++;
	}
      }
    }
    
    // Prefix sum to get the correct column pointers
    for (int i = 0; i < AAT_dim; ++i) 
      AAT.start_[i+1] += AAT.start_[i];
 
    AAT.index_.resize(AAT.start_.back());
    AAT.value_.resize(AAT.start_.back());
    AAT.p_end_ = AAT.start_;
    AAT.p_end_.back() = AAT.index_.size();
 
    std::vector<int> current_positions = AAT.start_;
 
    // Second pass to actually fill in the indices and values
    for (const auto& val : non_zero_values){
        int i = std::get<0>(val);
        int j = std::get<1>(val);
        double dot = std::get<2>(val);
 
        AAT.index_[current_positions[i]] = j;
        AAT.value_[current_positions[i]] = dot;
        current_positions[i]++;
        AAT.p_end_[i] = current_positions[i];
 
        if (i != j){
            AAT.index_[current_positions[j]] = i;
            AAT.value_[current_positions[j]] = dot;
            current_positions[j]++;
            AAT.p_end_[j] = current_positions[j];
        }
    }
    AAT.p_end_.clear();
    return AAT;
}

void newtonSolve(const HighsSparseMatrix &highs_a,
		 const std::vector<double> &theta,
                 const std::vector<double> &rhs,
		 std::vector<double> &lhs,
		 const int option_max_dense_col,
		 const double option_dense_col_tolerance) {
  HighsSparseMatrix AAT = computeAThetaAT_inner_product(highs_a, theta.data());

}
