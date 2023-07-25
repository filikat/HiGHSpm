#include "Direct.h"
#include "spral.h"
#include <cmath>

int newtonSolve(const HighsSparseMatrix &highs_a,
		const std::vector<double> &theta,
		const std::vector<double> &rhs,
		std::vector<double> &lhs,
		const int option_max_dense_col,
		const double option_dense_col_tolerance
		//		, ExperimentData& data
		) {

  /*
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

  HighsSparseMatrix AAT = computeAThetaAT(highs_a, theta.data());
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
  */
  return 1;
}
