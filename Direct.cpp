#include "Direct.h"
#include "spral.h"
#include <cmath>

double getWallTime() {
    using namespace std::chrono;
    using wall_clock = std::chrono::high_resolution_clock;
    return duration_cast<duration<double> >(wall_clock::now().time_since_epoch()).count();
}

double fillIn_LL(int nnz_AAT, int nnz_L, int MatrixSize) {
    return (double) (2*nnz_L - MatrixSize)/ nnz_AAT;
}

std::ostream& operator<<(std::ostream& os, const ExperimentData& data) {
    os << "decomposer: " <<data.decomposer << "\n model name:" << data.model_name << "\n model size: " << data.model_size << "\n nnz in AAT: " << data.nnz_AAT << "\n nnz in L: " << data.nnz_L << "\n solution error: " << data.solution_error << "\n residual error: " << data.residual_error << "\n fill-in: " << data.fill_in_factor << "\n analyse time: " << data.analysis_time<< "\n factorization time: " << data.factorization_time<< "\n solve time: " << data.solve_time << "\n time taken: " << data.time_taken << "\n";
    return os;
}

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

int newtonSolve(const HighsSparseMatrix &highs_a,
		const std::vector<double> &theta,
		const std::vector<double> &rhs,
		std::vector<double> &lhs,
		const int option_max_dense_col,
		const double option_dense_col_tolerance,
		ExperimentData& data) {

  int num_dense_col = 0;
  double max_density = 0; 
  for (int col = 0; col < highs_a.num_col_; col++) {
    double density = double(highs_a.start_[col+1] - highs_a.start_[col]) / double(highs_a.num_row_);
    if (density > option_dense_col_tolerance) {
      max_density = std::max(density, max_density);
      num_dense_col++;
    }
  }
  printf("Problem has %d dense columns (max = %g) for tolerance of %g\n",
	 int(num_dense_col), max_density, option_dense_col_tolerance);

  HighsSparseMatrix AAT = computeAThetaAT_inner_product(highs_a, theta.data());
  data.reset();
  data.decomposer = "ssids";
  data.model_size = highs_a.num_row_;
  data.nnz_AAT = AAT.numNz();
  // Prepare data structures for SPRAL
  std::vector<long> ptr;
  std::vector<int> row;
  std::vector<double> val;

  // Extract upper triangular part of AAT
  for (int col = 0; col < AAT.num_col_; col++){
    ptr.push_back(val.size());
    for (int idx = AAT.start_[col]; idx < AAT.start_[col+1]; idx++){
      int row_idx = AAT.index_[idx];
      if (row_idx >= col){
	val.push_back(AAT.value_[idx]);
	row.push_back(row_idx + 1);
      }
    }
  }

  for (auto& p : ptr) {
    ++p;
  }
  ptr.push_back(val.size() + 1);

  long* ptr_ptr = ptr.data();
  int* row_ptr = row.data();
  double* val_ptr = val.data();

  // Derived types
  void *akeep, *fkeep;
  struct spral_ssids_options options;
  struct spral_ssids_inform inform;

  // Initialize derived types
  akeep = NULL; fkeep = NULL;
  spral_ssids_default_options(&options);
  options.array_base = 1; // Need to set to 1 if using Fortran 1-based indexing 
  
  // Dense vector B 
  std::vector<double> B = rhs;

  // Perform analyse and factorise with data checking 
  bool check = true;
  double start_time = getWallTime();
  spral_ssids_analyse(check, AAT.num_col_, NULL, ptr_ptr, row_ptr, NULL, &akeep, &options, &inform);
  double end_time = getWallTime();
  data.analysis_time = end_time - start_time;
  if(inform.flag<0) {
    spral_ssids_free(&akeep, &fkeep);
    throw std::runtime_error("Error in spral_ssids_analyse");
  }
  
  bool positive_definite =true;
  double start_time2 = getWallTime();
  spral_ssids_factor(positive_definite, NULL, NULL, val_ptr, NULL, akeep, &fkeep, &options, &inform);
  double end_time2 = getWallTime();
  data.factorization_time = end_time2 - start_time2;
  if(inform.flag<0) {
    spral_ssids_free(&akeep, &fkeep);
    throw std::runtime_error("Error in spral_ssids_factor");
  }
  //Return the diagonal entries of the Cholesky factor
  std::vector<double> d(AAT.num_col_);
  
  void spral_ssids_enquire_posdef(const void *akeep,
				  const void *fkeep,
				  const struct spral_ssids_options *options,
				  struct spral_ssids_inform *inform,
				  double *d);
  if (inform.flag<0){
    spral_ssids_free(&akeep, &fkeep);
    throw std::runtime_error("Error in spral_ssids_enquire_posdef");
  }

  // Solve 
  double start_time3 = getWallTime();
  spral_ssids_solve1(0, B.data(), akeep, fkeep, &options, &inform);
  double end_time3 = getWallTime();
  data.solve_time = end_time3 - start_time3;
  if(inform.flag<0) {
    spral_ssids_free(&akeep, &fkeep);
    throw std::runtime_error("Error in spral_ssids_solve1");
  }
  data.nnz_L = inform.num_factor;
  data.time_taken = data.analysis_time + data.factorization_time + data.solve_time;
  //store solution in lhs
  lhs = B;
  // Free the memory allocated for SPRAL
  int cuda_error = spral_ssids_free(&akeep, &fkeep);
  if (cuda_error != 0){
    throw std::runtime_error("CUDA error in spral_ssids_free");
  }
  data.fill_in_factor = fillIn_LL(data.nnz_AAT, data.nnz_L, data.model_size);
  return 1;
}
