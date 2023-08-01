#include "Direct.h"
#include <cmath>
int IpmInvert::clear() {
  this->valid = false;
  this->system_size = -1;
  this->use_num_dense_col = 0;
  this->dense_col.clear();
  this->theta_d.clear();
  this->hatA.clear();
  return this->ssids_data.clear();
}

int SsidsData::clear() {
  // Free the memory allocated for SPRAL
  return spral_ssids_free(&this->akeep, &this->fkeep);
}

int augmentedInvert(const HighsSparseMatrix &highs_a,
		    const std::vector<double> &theta,
		    IpmInvert& invert,
		    ExperimentData& experiment_data) {
  assert(!invert.valid);
  
  assert(highs_a.isColwise());
  double start_time0 = getWallTime();
  experiment_data.reset();
  experiment_data.decomposer = "ssids";
  experiment_data.system_type = kSystemTypeAugmented;
  experiment_data.system_size = highs_a.num_col_ + highs_a.num_row_;
  experiment_data.system_nnz = highs_a.num_col_ + 2*highs_a.numNz();

  SsidsData& ssids_data = invert.ssids_data;
  int factor_status = callSsidsAugmentedFactor(highs_a, theta, ssids_data, experiment_data);
  experiment_data.time_taken = getWallTime() - start_time0;

  if (factor_status) return factor_status;
  invert.valid = true;
  return 0;
}

void augmentedSolve(const HighsSparseMatrix &highs_a,
		    const std::vector<double> &theta,
		    const std::vector<double> &rhs_x,
		    const std::vector<double> &rhs_y,
		    std::vector<double> &lhs_x,
		    std::vector<double> &lhs_y,
		    IpmInvert& invert,
		    ExperimentData& experiment_data) {
  assert(invert.valid);
  SsidsData& ssids_data = invert.ssids_data;
  
  double start_time = getWallTime();
  int row_index_offset = highs_a.num_col_;

  std::vector<double> rhs;
  for (int iCol = 0; iCol < highs_a.num_col_; iCol++)
    rhs.push_back(rhs_x[iCol]);
  for (int iRow = 0; iRow < highs_a.num_row_; iRow++)
    rhs.push_back(rhs_y[iRow]);
  
  callSsidsSolve(experiment_data.system_size, 1, rhs.data(), ssids_data);
  experiment_data.solve_time = getWallTime() - start_time;

  lhs_x.resize( highs_a.num_col_);
  lhs_y.resize( highs_a.num_row_);
  for (int iCol = 0; iCol < highs_a.num_col_; iCol++)
    lhs_x[iCol] = rhs[iCol];
  for (int iRow = 0; iRow < highs_a.num_row_; iRow++)
    lhs_y[iRow] = rhs[row_index_offset + iRow];

  experiment_data.residual_error = residualErrorAugmented(highs_a, theta, rhs_x, rhs_y, lhs_x, lhs_y);

}

int newtonInvert(const HighsSparseMatrix &highs_a,
		 const std::vector<double> &theta,
		 IpmInvert& invert,
		 const int option_max_dense_col,
		 const double option_dense_col_tolerance,
		 ExperimentData& experiment_data,
		 const bool quiet) {
  assert(!invert.valid);
  assert(highs_a.isColwise());
  const int system_size = highs_a.num_row_;
  std::vector<double> use_theta = theta;
  double use_dense_col_tolerance = option_dense_col_tolerance;
  //  use_dense_col_tolerance = 1.2;
  //  use_dense_col_tolerance = 0.1;

  double start_time0 = getWallTime();
  double start_time = start_time0;
  int col_max_nz = 0;
  std::vector<std::pair<double, int>> density_index;
  std::vector<double> analyse_density;

  std::vector<int>& dense_col = invert.dense_col;
  std::vector<double>& theta_d = invert.theta_d;
  std::vector<double>& hatA = invert.hatA;

  double max_sparse_col_density = 0;
  for (int iCol = 0; iCol < highs_a.num_col_; iCol++) {
    int col_nz = highs_a.start_[iCol+1] - highs_a.start_[iCol];
    double density_value = double(col_nz) / double(system_size);
    col_max_nz = std::max(col_nz, col_max_nz);
    analyse_density.push_back(density_value);
    if (density_value >= use_dense_col_tolerance) {
      density_index.push_back(std::make_pair(density_value, iCol));
    } else {
      max_sparse_col_density = std::max(density_value, max_sparse_col_density);
    }
  }
  const int model_num_dense_col = density_index.size();
  struct {
    bool operator()(std::pair<double, int> a, std::pair<double, int> b) const { return a.first > b.first; }
  }
  customMore;
  std::sort(density_index.begin(), density_index.end(), customMore);
  const int use_num_dense_col = std::min(model_num_dense_col, std::min(option_max_dense_col, system_size));
  
  // Take the first use_num_dense_col entries as dense
  for (int ix = 0; ix < use_num_dense_col; ix++) 
    dense_col.push_back(density_index[ix].second);
  if (use_num_dense_col < model_num_dense_col)
    max_sparse_col_density = density_index[use_num_dense_col].first;
  
  double max_density = double(col_max_nz)/ double(system_size);
  if (!quiet) {
    printf("Problem has %d rows and %d columns (max nonzeros = %d; density = %g) with %d dense at a tolerance of %g\n",
	   int(system_size), int(highs_a.num_col_),
	   int(col_max_nz), max_density, int(model_num_dense_col), use_dense_col_tolerance);
    analyseVectorValues(nullptr, "Column density", highs_a.num_col_, analyse_density);
  } else {
    //  printf("Newton solve uses %d dense columns\n", use_num_dense_col);
  }

  // Zero the entries of use_theta corresponding to dense columns
  for (int ix = 0; ix < use_num_dense_col; ix++) {
    int iCol = dense_col[ix];
    theta_d.push_back(theta[iCol]);
    use_theta[iCol] = 0;
  }

  HighsSparseMatrix AAT = computeAThetaAT(highs_a, use_theta);
  experiment_data.reset();
  experiment_data.decomposer = "ssids";
  experiment_data.system_type = kSystemTypeNewton;
  experiment_data.system_size = system_size;
  experiment_data.system_nnz = AAT.numNz();
  experiment_data.model_num_dense_col = model_num_dense_col;
  experiment_data.use_num_dense_col = use_num_dense_col;
  experiment_data.dense_col_tolerance = use_dense_col_tolerance;
  if (model_num_dense_col) assert(density_index[0].first == max_density);
  experiment_data.model_max_dense_col = max_density;
  experiment_data.system_max_dense_col = max_sparse_col_density;

  experiment_data.form_time = getWallTime() - start_time;

  SsidsData& ssids_data = invert.ssids_data;
  int factor_status = callSsidsNewtonFactor(AAT, ssids_data, experiment_data);
  if (factor_status) return factor_status;

  if (use_num_dense_col) {
    double start_time = getWallTime();
    // Set up RHS for use_num_dense_col columns
    hatA.assign(use_num_dense_col*system_size, 0);
    // First form \hat{A}_d for the dense columns
    int offset = 0;
    for (int ix = 0; ix < use_num_dense_col; ix++) {
      int iCol = dense_col[ix];
      for (int iEl = highs_a.start_[iCol]; iEl < highs_a.start_[iCol+1]; iEl++) 
	hatA[offset+highs_a.index_[iEl]] = highs_a.value_[iEl];
      offset += system_size;
    }
    callSsidsSolve(system_size, use_num_dense_col, hatA.data(), ssids_data);
    experiment_data.factorization_time += getWallTime() - start_time;
  }
  invert.system_size = system_size;
  invert.use_num_dense_col = use_num_dense_col;
  experiment_data.time_taken = getWallTime() - start_time0;
  invert.valid = true;
  return 0;
}

int newtonSolve(const HighsSparseMatrix &highs_a,
		const std::vector<double> &theta,
		const std::vector<double> &rhs,
		std::vector<double> &lhs,
		IpmInvert& invert,
		ExperimentData& experiment_data) {
  assert(invert.valid);
  SsidsData& ssids_data = invert.ssids_data;
  
  double start_time = getWallTime();

  const int system_size = invert.system_size;
  const int use_num_dense_col = invert.use_num_dense_col;
  const std::vector<int>& dense_col = invert.dense_col;
  const std::vector<double>& theta_d = invert.theta_d;
  const std::vector<double>& hatA = invert.hatA;
  lhs = rhs;
  callSsidsSolve(system_size, 1, lhs.data(), ssids_data);
  if (use_num_dense_col) {
    // Now form D = \Theta_d^{-1} + \hat{A}_d^TA_d and \hat_b = \hat{A}_d^Tb
    
    std::vector<std::vector<double>> d_matrix;
    std::vector<double> d_rhs;
    std::vector<double> d_sol;
    d_matrix.resize(use_num_dense_col);
    int offset = 0;
    for (int d_col = 0; d_col < use_num_dense_col; d_col++) {
      d_matrix[d_col].resize(use_num_dense_col);
      for (int d_row = 0; d_row < use_num_dense_col; d_row++) {
	int iCol = dense_col[d_row];
	double value = 0;
	for (int iEl = highs_a.start_[iCol]; iEl < highs_a.start_[iCol+1]; iEl++)
	  value += hatA[offset+highs_a.index_[iEl]] * highs_a.value_[iEl];
	d_matrix[d_col][d_row] = value;
      }
      d_matrix[d_col][d_col] += 1/theta_d[d_col];
      double value = 0;
      for (int iRow = 0; iRow < system_size; iRow++) 
	value += hatA[offset+iRow] * rhs[iRow];
      d_rhs.push_back(value);
      offset += system_size;
    }
    int gepp_status = gepp(d_matrix, d_rhs, d_sol);
    if (gepp_status) return 1;
    
    // Subtract \hat{A}_d^Td_sol from lhs
    offset = 0;
    for (int d_col = 0; d_col < use_num_dense_col; d_col++) {
      for (int iRow = 0; iRow < system_size; iRow++) 
	lhs[iRow] -= hatA[offset+iRow] * d_sol[d_col];
      offset += system_size;
    }
  }
  experiment_data.solve_time = getWallTime() - start_time;

  experiment_data.residual_error = residualErrorNewton(highs_a, theta, rhs, lhs);
  return 0;
}

bool increasingIndex(const HighsSparseMatrix& matrix) {
  if (matrix.isRowwise()) {
    for (int iRow = 0; iRow < matrix.num_row_; iRow++)
      for (int iEl = matrix.start_[iRow]+1; iEl < matrix.start_[iRow+1]; iEl++) 
	if (matrix.index_[iEl] <= matrix.index_[iEl-1]) return false;
  } else {
    for (int iCol = 0; iCol < matrix.num_col_; iCol++)
      for (int iEl = matrix.start_[iCol]+1; iEl < matrix.start_[iCol+1]; iEl++)
	if (matrix.index_[iEl] <= matrix.index_[iEl-1]) return false;
  }
  return true;
}

void productAThetaAT(const HighsSparseMatrix& matrix,
		     const std::vector<double>& theta,
		     const std::vector<double>& x,
		     std::vector<double>& result) {
  assert(int(x.size()) == matrix.num_row_);
  std::vector<double> ATx;
  matrix.productTranspose(ATx, x);
  if (!theta.empty()) 
    for (int ix = 0; ix < matrix.num_col_; ix++) ATx[ix] *= theta[ix];
  matrix.product(result, ATx);
}

HighsSparseMatrix computeAThetaAT(const HighsSparseMatrix& matrix,
				  const std::vector<double>& theta) {
  const bool scatter = true;
  // Create a row-wise copy of the matrix
  HighsSparseMatrix AT = matrix;
  AT.ensureRowwise();

  int AAT_dim = matrix.num_row_;
  HighsSparseMatrix AAT;
  AAT.num_col_ = AAT_dim;
  AAT.num_row_ = AAT_dim;
  AAT.start_.resize(AAT_dim+1,0);

  std::vector<std::tuple<int, int, double>> non_zero_values;

  // First pass to calculate the number of non-zero elements in each column
  //
  if (scatter) {
    std::vector<double> matrix_row(matrix.num_col_, 0);
    for (int iRow = 0; iRow < AAT_dim; iRow++) {
      for (int iEl = AT.start_[iRow]; iEl < AT.start_[iRow+1]; iEl++) 
	matrix_row[AT.index_[iEl]] = AT.value_[iEl];
      for (int iCol = iRow; iCol < AAT_dim; iCol++) {
	double dot = 0.0;
	for (int iEl = AT.start_[iCol]; iEl < AT.start_[iCol+1]; iEl++) {
	  const int ix = AT.index_[iEl];
	  const double theta_i = !theta.empty() ? theta[ix] : 1;
	  dot += theta_i * matrix_row[ix] * AT.value_[iEl];
	}
	if (dot != 0.0) {
	  non_zero_values.emplace_back(iRow, iCol, dot);
	  AAT.start_[iRow+1]++;
	  if (iRow != iCol) AAT.start_[iCol+1]++;
	}
      }
      for (int iEl = AT.start_[iRow]; iEl < AT.start_[iRow+1]; iEl++) 
	matrix_row[AT.index_[iEl]] = 0;
      for (int ix = 0; ix < matrix.num_col_; ix++)
	assert(!matrix_row[ix]);
    }
  } else {
    assert(increasingIndex(AT));
    for (int i = 0; i < AAT_dim; ++i) {
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
	    const double theta_i = !theta.empty() ? theta[AT.index_[k]] : 1;
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

// Gaussian elimination with partial pivoting for a dense matrix
//
// Returns 0 if it's successful
//
// Returns -1 if GE fails due to pivot less than kMinAbsPivot

const double kMinAbsPivot = 1e-12;
int gepp(const std::vector<std::vector<double>>& matrix,
	 const std::vector<double>& rhs,
	 std::vector<double>& solution) {
  std::vector<std::vector<double>> ge_matrix = matrix;
  solution = rhs;
  int dim = int(rhs.size());
  for (int k = 0; k < dim; k++) {
    double max_pivot_col_value = kMinAbsPivot;
    int pivot_row = -1;
    for (int l = k; l < dim; l++) {
      double abs_pivot_col_value = std::abs(ge_matrix[l][k]);
      if (max_pivot_col_value < abs_pivot_col_value) {
	max_pivot_col_value = abs_pivot_col_value;
	pivot_row = l;
      }
    }
    // No pivot so matrix is singular
    if (pivot_row < 0) return -1;
    if (pivot_row != k) {
      // Perform row interchange
      for (int iCol = k; iCol < dim; iCol++) {
	double row_k_v = ge_matrix[k][iCol];
	ge_matrix[k][iCol] = ge_matrix[pivot_row][iCol];
	ge_matrix[pivot_row][iCol] = row_k_v;
      }
      double rhs_k_v = solution[k];
      solution[k] = solution[pivot_row];
      solution[pivot_row] = rhs_k_v;
    }
    assert(ge_matrix[k][k] != 0);
    for (int iRow = k+1; iRow < dim; iRow++) {
      double multiplier = ge_matrix[iRow][k] / ge_matrix[k][k] ;
      assert(std::fabs(multiplier) <= 1);
      for (int iCol = k+1; iCol < dim; iCol++) 
	ge_matrix[iRow][iCol] -= multiplier * ge_matrix[k][iCol];
      solution[iRow] -= multiplier * solution[k];
    }
  }
  for (int iRow = dim-1; iRow >= 0; iRow--) {
    for (int iCol = iRow+1; iCol < dim; iCol++)
      solution[iRow] -= solution[iCol] * ge_matrix[iRow][iCol];
    solution[iRow] /= ge_matrix[iRow][iRow];
  }
  return 0;
}

//int ssids_decompose();

int callSsidsAugmentedFactor(const HighsSparseMatrix& matrix,
			     const std::vector<double>& theta,
			     SsidsData& ssids_data,
			     ExperimentData& experiment_data) {
  experiment_data.form_time = 0;
  double start_time = getWallTime();

  // Prepare data structures for SPRAL
  const int array_base = 0;
  std::vector<long> ptr;
  std::vector<int> row;
  std::vector<double> val;
  spral_ssids_default_options(&ssids_data.options);

  // Extract lower triangular part of augmented system
  //
  // First the columns for [-\Theta^{-1}]
  //                       [      A     ]
  int row_index_offset = matrix.num_col_;
  for (int iCol = 0; iCol < matrix.num_col_; iCol++) {
    const double theta_i = !theta.empty() ? theta[iCol] : 1;
    ptr.push_back(val.size());
    val.push_back(-1/theta_i);
    row.push_back(iCol + array_base);
    for (int iEl = matrix.start_[iCol]; iEl < matrix.start_[iCol+1]; iEl++) {
      int iRow = matrix.index_[iEl];
      val.push_back(matrix.value_[iEl]);
      row.push_back(row_index_offset + iRow + array_base);
    }
  }
  // Now the (zero) columns for [A^T]
  //                            [ 0 ]
  
  const double diagonal = ssids_data.options.small;//1e-20;
  for (int iRow = 0; iRow < matrix.num_row_; iRow++) {
    ptr.push_back(val.size());
    if (diagonal) {
      val.push_back(diagonal);
      row.push_back(row_index_offset + iRow + array_base);
    }
  }
  // Possibly add 1 to all the pointers
  if (array_base) for (auto& p : ptr) ++p;

  // Add the last pointer
  ptr.push_back(val.size() + array_base);
 
  long* ptr_ptr = ptr.data();
  int* row_ptr = row.data();
  double* val_ptr = val.data();

  // Initialize derived types
  ssids_data.akeep = nullptr;
  ssids_data.fkeep = nullptr;
  // Need to set to 1 if using Fortran 1-based indexing 
  ssids_data.options.array_base = array_base; 
  //  ssids_data.options.print_level = 2; 
  
  experiment_data.setup_time = getWallTime() - start_time;

  // Perform analyse and factorise with data checking 
  bool check = true;
  start_time = getWallTime();
  spral_ssids_analyse(check, experiment_data.system_size,
		      nullptr, ptr_ptr, row_ptr, nullptr,
		      &ssids_data.akeep,
		      &ssids_data.options,
		      &ssids_data.inform);
  experiment_data.analysis_time = getWallTime() - start_time;
  if (ssids_data.inform.flag < 0) return 1;
  
  bool positive_definite = false;
  start_time = getWallTime();
  spral_ssids_factor(positive_definite,
		     nullptr, nullptr, val_ptr, nullptr,
		     ssids_data.akeep,
		     &ssids_data.fkeep,
		     &ssids_data.options,
		     &ssids_data.inform);
  experiment_data.factorization_time = getWallTime() - start_time;
  if(ssids_data.inform.flag < 0) return 1;
  experiment_data.nnz_L = ssids_data.inform.num_factor;
  experiment_data.fillIn_LL();
  return 0;
}

int callSsidsNewtonFactor(const HighsSparseMatrix& AThetaAT,
			  SsidsData& ssids_data,
			  ExperimentData& experiment_data) {
  double start_time = getWallTime();
  // Prepare data structures for SPRAL
  const int array_base = 0;
  std::vector<long> ptr;
  std::vector<int> row;
  std::vector<double> val;
  spral_ssids_default_options(&ssids_data.options);

  // Extract lower triangular part of AAT
  for (int col = 0; col < AThetaAT.num_col_; col++){
    ptr.push_back(val.size());
    for (int idx = AThetaAT.start_[col]; idx < AThetaAT.start_[col+1]; idx++){
      int row_idx = AThetaAT.index_[idx];
      if (row_idx >= col){
	val.push_back(AThetaAT.value_[idx]);
	row.push_back(row_idx + array_base);
      }
    }
  }
  // Possibly add 1 to all the pointers
  if (array_base) for (auto& p : ptr) ++p;

  // Add the last pointer
  ptr.push_back(val.size() + array_base);

  long* ptr_ptr = ptr.data();
  int* row_ptr = row.data();
  double* val_ptr = val.data();

  // Initialize derived types
  ssids_data.akeep = nullptr;
  ssids_data.fkeep = nullptr;
  // Need to set to 1 if using Fortran 1-based indexing 
  ssids_data.options.array_base = array_base; // Need to set to 1 if using Fortran 1-based indexing 

  experiment_data.setup_time = getWallTime() - start_time;
  
  // Perform analyse and factorise with data checking 
  bool check = true;
  start_time = getWallTime();
  spral_ssids_analyse(check, AThetaAT.num_col_,
		      nullptr, ptr_ptr, row_ptr, nullptr,
		      &ssids_data.akeep,
		      &ssids_data.options,
		      &ssids_data.inform);
  experiment_data.analysis_time = getWallTime() - start_time;
  if(ssids_data.inform.flag < 0) return 1;
  
  bool positive_definite =true;
  start_time = getWallTime();
  spral_ssids_factor(positive_definite,
		     nullptr, nullptr, val_ptr, nullptr,
		     ssids_data.akeep,
		     &ssids_data.fkeep,
		     &ssids_data.options,
		     &ssids_data.inform);
  experiment_data.factorization_time = getWallTime() - start_time;
  if(ssids_data.inform.flag < 0) return 1;

  /*
  //Return the diagonal entries of the Cholesky factor
  std::vector<double> d(AThetaAT.num_col_);
  void spral_ssids_enquire_posdef(const void *akeep,
				  const void *fkeep,
				  const struct spral_ssids_options *options,
				  struct spral_ssids_inform *inform,
				  double *d);
  */
  experiment_data.nnz_L = ssids_data.inform.num_factor;
  experiment_data.fillIn_LL();
  return 0;
}

void callSsidsSolve(const int system_size,
		    const int num_rhs,
		    double* rhs,
		    SsidsData& ssids_data) {
  spral_ssids_solve(0, num_rhs, rhs,
		    system_size,
		    ssids_data.akeep,
		    ssids_data.fkeep,
		    &ssids_data.options,
		    &ssids_data.inform);
}
