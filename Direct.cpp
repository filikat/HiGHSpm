#include "Direct.h"
#include <cmath>
void IpmInvert::clear() {
  this->valid = false;
  this->system_size = -1;
  this->use_num_dense_col = 0;
  this->dense_col.clear();
  this->theta_d.clear();
  this->hatA_d.clear();
  this->d_matrix.clear();
  const int decomposer_source = this->decomposer_source;
  if (decomposer_source == kDecomposerSourceSsids){
    this->ssids_data.clear();
  } else if (decomposer_source == kDecomposerSourceMa86){
    this->ma86_data.clear();
  } else if (decomposer_source == kDecomposerSourceQdldl){
    this->qdldl_data.clear();
  } else if (decomposer_source == kDecomposerSourceCholmod){
    this->cholmod_data.clear();
  } else if (decomposer_source == kDecomposerSourceHighs){
    //    this->highs_factor.clear();
  }
}

int SsidsData::clear() {
  // Free the memory allocated for SPRAL
#ifdef HAVE_SPRAL
  if (spral_ssids_free(&this->akeep, &this->fkeep)) return kDecomposerStatusErrorClear;
  return kDecomposerStatusOk;
#else
  return kDecomposerStatusErrorClear;
#endif
}

void CholmodData::clear() {   
#ifdef HAVE_CHOLMOD
  cholmod_free_factor(&this->L, &this->c);
  cholmod_free_sparse(&this->a, &this->c);
  cholmod_free_triplet(&this->T, &this->c);
  //cholmod_free_dense (&r, &c) ;
  cholmod_free_dense(&this->x, &this->c);
  cholmod_free_dense(&this->b, &this->c);
  //cholmod_free_sparse(&L_sparse, &c);
  cholmod_finish(&c);  
#endif
}

void chooseDenseColumns(const HighsSparseMatrix &highs_a,
			const std::vector<double> &theta, 
			const int option_max_dense_col,
			const double option_dense_col_tolerance,
			std::vector<int> &dense_col,
			ExperimentData &experiment_data,
			const bool quiet) {
  double use_dense_col_tolerance = option_dense_col_tolerance;
  //  use_dense_col_tolerance = 1.2;
  //  use_dense_col_tolerance = 0.1;
  const int system_size = experiment_data.system_size;
  int col_max_nz = 0;
  std::vector<std::pair<double, int>> density_index;
  std::vector<double> analyse_density;

  double max_sparse_col_density = 0;
  for (int iCol = 0; iCol < highs_a.num_col_; iCol++) {
    int col_nz = highs_a.start_[iCol + 1] - highs_a.start_[iCol];
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
    bool operator()(std::pair<double, int> a, std::pair<double, int> b) const {
      return a.first > b.first;
    }
  } customMore;
  std::sort(density_index.begin(), density_index.end(), customMore);
  // Number of dense columns to be used cannot exceed the number of dense columns permitted...
  int use_num_dense_col = option_max_dense_col;
  // ... the number of dense columns in the model...
  use_num_dense_col = std::min(model_num_dense_col, use_num_dense_col);
  // ... and must leave at least as many sparse columns as the system
  // size
  use_num_dense_col = std::min(highs_a.num_col_-system_size, use_num_dense_col);
  if (use_num_dense_col < model_num_dense_col)
    max_sparse_col_density = density_index[use_num_dense_col].first;

  std::vector<bool> is_dense(highs_a.num_col_, false);
  // Take the first use_num_dense_col entries as dense, counting how many 
  for (int ix = 0; ix < use_num_dense_col; ix++) {
    int iCol = density_index[ix].second;
    dense_col.push_back(iCol);
    is_dense[iCol] = true;
  }
  // Find the largest theta value amongst the sparse columns
  double max_sparse_theta = 0;
  for (int iCol = 0; iCol < highs_a.num_col_; iCol++) {
    if (is_dense[iCol]) continue;
    max_sparse_theta = std::max(theta[iCol], max_sparse_theta);
  }
  const bool check_num_ok_sparse_theta = false;//true;//
  if (check_num_ok_sparse_theta) {
    double ok_sparse_theta = max_sparse_theta * ok_theta_relative_tolerance;
    for (;;) {
      // Count the number of sparse theta values that are at least ok_sparse_theta
      int num_ok_sparse_theta = 0;
      for (int iCol = 0; iCol < highs_a.num_col_; iCol++) {
	if (is_dense[iCol]) continue;
	if (theta[iCol] > ok_sparse_theta) num_ok_sparse_theta++;
      }
      if (num_ok_sparse_theta >= system_size) break;
      // Danger of bad numerics since there are fewer OK theta values in
      // sparse columns than the number of rows
      std::vector<int>new_dense_col;
      // Work through the dense columns from low to high density,
      // removing those with OK theta values until there are sufficient
      // sparse theta values that are at least ok_sparse_theta
      for (int ix = use_num_dense_col-1; ix >=0 ; ix--) {
	int iCol = dense_col[ix];
	if (theta[iCol] > ok_sparse_theta &&
	    num_ok_sparse_theta < system_size) {
	  // Add this column to the sparse columns
	  is_dense[iCol] = false;
	  num_ok_sparse_theta++;
	  // Update the record of the densest sparse column
	  int col_nz = highs_a.start_[iCol + 1] - highs_a.start_[iCol];
	  double density_value = double(col_nz) / double(system_size);
	  max_sparse_col_density = std::max(density_value, max_sparse_col_density);
	  printf("Eliminate dense column %d of density %g\n", iCol, density_value);
	} else {
	  new_dense_col.push_back(iCol);
	}
      }
      dense_col = new_dense_col;
      use_num_dense_col = dense_col.size();
      if (num_ok_sparse_theta >= system_size) break;
      // Still not enough ok_sparse_theta, so repeat with smaller
      // requirement for OK theta
      ok_sparse_theta /= 10;
    }
    const bool check_method = true;
    if (check_method) {
      // Check that algorithm has worked
      int num_ok_sparse_theta = 0;
      for (int iCol = 0; iCol < highs_a.num_col_; iCol++) {
	if (is_dense[iCol]) continue;
	if (theta[iCol] > ok_sparse_theta) num_ok_sparse_theta++;
      }
      assert(num_ok_sparse_theta >= system_size);
    }
  }
  
  double max_density = double(col_max_nz) / double(system_size);
  if (!quiet) {
    printf("Problem has %d rows and %d columns (max nonzeros = %d; density = "
           "%g) with %d dense at a tolerance of %g\n",
           int(system_size), int(highs_a.num_col_), int(col_max_nz),
           max_density, int(model_num_dense_col), use_dense_col_tolerance);
    analyseVectorValues(nullptr, "Column density", highs_a.num_col_,
                        analyse_density);
  }

  experiment_data.model_num_dense_col = model_num_dense_col;
  experiment_data.use_num_dense_col = use_num_dense_col;
  experiment_data.dense_col_tolerance = use_dense_col_tolerance;
  if (model_num_dense_col)
    assert(density_index[0].first == max_density);
  experiment_data.model_max_dense_col = max_density;
  experiment_data.system_max_dense_col = max_sparse_col_density;

}

void MA86Data::clear() {
  // Free the memory allocated for MA86
#ifdef HAVE_MA86
  wrapper_ma86_finalise(&this->keep, &this->control);
#endif
}

void QDLDLData::clear() {
#ifdef HAVE_QDLDL
  free(this->Lp);
  free(this->Li);
  free(this->Lx);
  free(this->D);
  free(this->Dinv);
  free(this->etree);
  free(this->Lnz);
  free(this->iwork);
  free(this->bwork);
  free(this->fwork);
  free(this->x);
#endif
}

/*
1-Spral
2-MA86
3-QDLQL
5-HiGHS
*/
int augmentedInvert(const HighsSparseMatrix &highs_a,
                    const std::vector<double> &theta, IpmInvert &invert,
                    ExperimentData &experiment_data,
		    const int decomposer_source) {
  assert(!invert.valid);
  assert(highs_a.isColwise());

  invert.decomposer_source = decomposer_source;

  double start_time0 = getWallTime();
  assert(decomposer_source >= kDecomposerSourceMin &&
	 decomposer_source <= kDecomposerSourceMax &&
	 decomposer_source != kDecomposerSourceCholmod);

  int factor_status;
  experiment_data.system_type = kSystemTypeAugmented;
  experiment_data.system_size = highs_a.num_col_ + highs_a.num_row_;
  if (decomposer_source == kDecomposerSourceSsids) {
    experiment_data.decomposer = "ssids";
    factor_status =
      callSsidsAugmentedFactor(highs_a, theta, invert.ssids_data, experiment_data);
  } else if (decomposer_source == kDecomposerSourceMa86) {
    experiment_data.decomposer = "ma86";
    factor_status = 
      callMA86AugmentedFactor(highs_a, theta, invert.ma86_data, experiment_data);
  } else if (decomposer_source == kDecomposerSourceQdldl){
    experiment_data.decomposer = "qdldl";
    factor_status = 
      callQDLDLAugmentedFactor(highs_a, theta, invert.qdldl_data, experiment_data);
  } else if (decomposer_source == kDecomposerSourceHighs) {
    experiment_data.decomposer = "qdldl";
    factor_status = 
      callHighsAugmentedFactor(highs_a, theta, invert.highs_factor, experiment_data);
  }
  experiment_data.nla_time.total = getWallTime() - start_time0;
  if (factor_status) return factor_status;
  invert.valid = true;
  return kDecomposerStatusOk;    
}

void augmentedSolve(const HighsSparseMatrix &highs_a,
                    const std::vector<double> &theta,
                    const std::vector<double> &rhs_x,
                    const std::vector<double> &rhs_y,
                    std::vector<double> &lhs_x, std::vector<double> &lhs_y,
                    IpmInvert &invert, ExperimentData &experiment_data,
                    const int decomposer_source) {
  assert(invert.valid);
  SsidsData &ssids_data = invert.ssids_data;
  MA86Data &ma86_data = invert.ma86_data;
  QDLDLData &qdldl_data = invert.qdldl_data;

  double start_time = getWallTime();
  int row_index_offset = highs_a.num_col_;

  std::vector<double> rhs;
  for (int iCol = 0; iCol < highs_a.num_col_; iCol++)
    rhs.push_back(rhs_x[iCol]);
  for (int iRow = 0; iRow < highs_a.num_row_; iRow++)
    rhs.push_back(rhs_y[iRow]);
  int system_size = highs_a.num_col_ + highs_a.num_row_;

  assert(invert.decomposer_source == decomposer_source);
  if (decomposer_source == kDecomposerSourceSsids) {
    callSsidsSolve(system_size, 1, rhs.data(), ssids_data);
  } else if (decomposer_source == kDecomposerSourceMa86) {
    callMA86Solve(system_size, 1, rhs.data(), ma86_data);
  } else if (decomposer_source == kDecomposerSourceQdldl) {
    callQDLDLSolve(system_size, 1, rhs.data(), invert.qdldl_data);
  } else if (decomposer_source == kDecomposerSourceHighs) {
    callHighsSolve(system_size, 1, rhs.data(), invert.highs_factor);
  }
  experiment_data.nla_time.solve = getWallTime() - start_time;

  lhs_x.resize(highs_a.num_col_);
  lhs_y.resize(highs_a.num_row_);
  for (int iCol = 0; iCol < highs_a.num_col_; iCol++)
    lhs_x[iCol] = rhs[iCol];
  for (int iRow = 0; iRow < highs_a.num_row_; iRow++)
    lhs_y[iRow] = rhs[row_index_offset + iRow];

  experiment_data.residual_error =
      residualErrorAugmented(highs_a, theta, rhs_x, rhs_y, lhs_x, lhs_y);
}

double augmentedCondition(const HighsSparseMatrix &matrix,
                          const std::vector<double> &theta, IpmInvert &invert,
                          const int decomposer_source) {
  assert(invert.valid);
  SsidsData &ssids_data = invert.ssids_data;
  const bool chatty = false;

  // Get the largest eigenvalue by power method
  const double x_dim = matrix.num_col_;
  const double y_dim = matrix.num_row_;
  std::vector<double> from_x_iterate(x_dim);
  std::vector<double> to_x_iterate(x_dim);
  std::vector<double> from_y_iterate(y_dim);
  std::vector<double> to_y_iterate(y_dim);
  HighsRandom random;
  double lambda_n = 0;
  for (int iRow = 0; iRow < x_dim; iRow++)
    from_x_iterate[iRow] = random.fraction();
  for (int iRow = 0; iRow < y_dim; iRow++)
    from_y_iterate[iRow] = random.fraction();
  for (int iter = 0; iter < 10; iter++) {
    double from_iterate_norm = Norm2(from_x_iterate, from_y_iterate);
    assert(from_iterate_norm > 0);
    VectorScale(from_x_iterate, 1 / from_iterate_norm);
    VectorScale(from_y_iterate, 1 / from_iterate_norm);

    std::vector<double> ATy;
    matrix.productTranspose(ATy, from_y_iterate);
    std::vector<double> Ax;
    matrix.product(Ax, from_x_iterate);
    for (int ix = 0; ix < x_dim; ix++) {
      const double theta_i = !theta.empty() ? theta[ix] : 1;
      to_x_iterate[ix] = -1 / theta_i * from_x_iterate[ix] + ATy[ix];
    }
    to_y_iterate = Ax;

    double to_iterate_norm = Norm2(to_x_iterate, to_y_iterate);
    const bool converged =
        std::fabs(lambda_n - to_iterate_norm) / std::max(1.0, lambda_n) <
            1e-2 &&
        iter > 2;
    lambda_n = to_iterate_norm;
    if (chatty)
      printf("Iter %2d: lambda_n = %11.4g\n", iter, lambda_n);
    if (converged)
      break;
    from_x_iterate = to_x_iterate;
    from_y_iterate = to_y_iterate;
  }

  ExperimentData experiment_data;
  double lambda_1 = 0;
  for (int iRow = 0; iRow < x_dim; iRow++)
    from_x_iterate[iRow] = random.fraction();
  for (int iRow = 0; iRow < y_dim; iRow++)
    from_y_iterate[iRow] = random.fraction();
  for (int iter = 0; iter < 10; iter++) {
    double from_iterate_norm = Norm2(from_x_iterate, from_y_iterate);
    assert(from_iterate_norm > 0);
    VectorScale(from_x_iterate, 1 / from_iterate_norm);
    VectorScale(from_y_iterate, 1 / from_iterate_norm);

    augmentedSolve(matrix, theta, from_x_iterate, from_y_iterate, to_x_iterate,
                   to_y_iterate, invert, experiment_data, decomposer_source);

    double to_iterate_norm = Norm2(to_x_iterate, to_y_iterate);
    const bool converged =
        std::fabs(lambda_1 - to_iterate_norm) / std::max(1.0, lambda_1) <
            1e-2 &&
        iter > 2;
    lambda_1 = to_iterate_norm;
    if (chatty)
      printf("Iter %2d: lambda_1 = %11.4g\n", iter, lambda_1);
    if (converged)
      break;
    from_x_iterate = to_x_iterate;
    from_y_iterate = to_y_iterate;
  }
  double condition = lambda_n / lambda_1;
  if (chatty)
    printf("Condition = %g\n", condition);

  return condition;
}

int newtonInvert(const HighsSparseMatrix &highs_a,
                 const std::vector<double> &theta, IpmInvert &invert,
                 const int option_max_dense_col,
                 const double option_dense_col_tolerance,
                 ExperimentData &experiment_data, const bool quiet,
		             const int decomposer_source) {
  const bool first_call_with_theta = !invert.valid;
  assert(first_call_with_theta);
  assert(highs_a.isColwise());
  const int system_size = highs_a.num_row_;
  std::vector<double> use_theta = theta;

  double start_time0 = getWallTime();
  double start_time = start_time0;

  std::vector<int> &dense_col = invert.dense_col;
  std::vector<double> &theta_d = invert.theta_d;
  std::vector<double> &hatA_d = invert.hatA_d;

  experiment_data.system_type = kSystemTypeNewton;
  experiment_data.system_size = system_size;

  chooseDenseColumns(highs_a, theta,
		     option_max_dense_col, option_dense_col_tolerance,
		     dense_col, experiment_data, quiet);
  int use_num_dense_col = dense_col.size();
  // Zero the entries of use_theta corresponding to dense columns  
  for (int ix = 0; ix < use_num_dense_col; ix++) {
    int iCol = dense_col[ix];
    theta_d.push_back(theta[iCol]);
    use_theta[iCol] = 0;
  }
  // Possibly zero the entries of use_theta corresponding to small
  // values of theta
  const double zero_theta_relative_tolerance = 0;//1e-8;//
  const double zero_theta_tolerance =  experiment_data.theta_max * zero_theta_relative_tolerance;
  if (zero_theta_tolerance>0) {
    int num_zeroed_theta = 0;
    for (int iCol = 0; iCol < highs_a.num_col_; iCol++) {
      if (std::fabs(use_theta[iCol]) < zero_theta_tolerance) {
	use_theta[iCol] = 0;
	num_zeroed_theta++;
      }
    }
    if (num_zeroed_theta) printf("Zeroed %d/%d theta values less than %g\n",
				 num_zeroed_theta,
				 highs_a.num_col_,
				 zero_theta_tolerance);
  }
  HighsSparseMatrix AAT;
  int AAT_status = computeAThetaAT(highs_a, use_theta, AAT);
  if (AAT_status) return AAT_status;

  const bool aat_quiet = false;
  if (!aat_quiet) 
    analyseVectorValues(nullptr, "A.Theta.A^T nonzeros", AAT.numNz(),
			AAT.value_);

  experiment_data.system_type = kSystemTypeNewton;
  experiment_data.system_size = system_size;
  experiment_data.system_nnz = AAT.numNz();

  experiment_data.nla_time.form = getWallTime() - start_time;
  SsidsData &ssids_data = invert.ssids_data;
  MA86Data &ma86_data = invert.ma86_data;
  QDLDLData &qdldl_data = invert.qdldl_data;
  CholmodData &cholmod_data = invert.cholmod_data;
  int factor_status;
  assert(decomposer_source >= kDecomposerSourceMin &&
	 decomposer_source <= kDecomposerSourceMax);
  if (decomposer_source == kDecomposerSourceSsids){
    experiment_data.decomposer = "ssids";
    factor_status = callSsidsNewtonFactor(AAT, ssids_data, experiment_data);
  } else if (decomposer_source == kDecomposerSourceMa86) {
    experiment_data.decomposer = "ma86";
    factor_status = callMA86NewtonFactor(AAT, ma86_data, experiment_data);
  } else if (decomposer_source == kDecomposerSourceQdldl) {
    experiment_data.decomposer = "qdldl";
    factor_status = callQDLDLNewtonFactor(AAT, qdldl_data, experiment_data);
  } else if (decomposer_source == kDecomposerSourceCholmod) {
    experiment_data.decomposer = "cholmod";
    factor_status = callCholmodNewtonFactor(AAT, cholmod_data, experiment_data);
  } else if (decomposer_source == kDecomposerSourceHighs) {
    experiment_data.decomposer = "HiGHS";
    assert(111==222);
    //    factor_status = callHighsNewtonFactor(AAT, Highs_data, experiment_data);
  }

  if (factor_status)
    return factor_status;

  if (use_num_dense_col) {
    double start_time = getWallTime();
    // Set up RHS for use_num_dense_col columns
    hatA_d.assign(use_num_dense_col * system_size, 0);
    // First form \hat{A}_d for the dense columns
    int offset = 0;
    for (int ix = 0; ix < use_num_dense_col; ix++) {
      int iCol = dense_col[ix];
      for (int iEl = highs_a.start_[iCol]; iEl < highs_a.start_[iCol + 1];
           iEl++)
        hatA_d[offset + highs_a.index_[iEl]] = highs_a.value_[iEl];
      offset += system_size;
    }
    if (decomposer_source == kDecomposerSourceSsids){
      callSsidsSolve(system_size, use_num_dense_col, hatA_d.data(), ssids_data);
    } else if (decomposer_source == kDecomposerSourceMa86) {
      callMA86Solve(system_size, use_num_dense_col, hatA_d.data(), ma86_data);
    } else if (decomposer_source == kDecomposerSourceQdldl){
      callQDLDLSolve(system_size, use_num_dense_col, hatA_d.data(), qdldl_data);
    } else if (decomposer_source == kDecomposerSourceCholmod) {
      callCholmodSolve(system_size, use_num_dense_col, hatA_d.data(), cholmod_data);
    } else if (decomposer_source == kDecomposerSourceHighs){
      callHighsSolve(system_size, use_num_dense_col, hatA_d.data(), invert.highs_factor);
    }
    
    // Now form D = \Theta_d^{-1} + \hat{A}_d^TA_d
    std::vector<std::vector<double>> &d_matrix = invert.d_matrix;
    d_matrix.resize(use_num_dense_col);
    offset = 0;
    for (int d_col = 0; d_col < use_num_dense_col; d_col++) {
      d_matrix[d_col].resize(use_num_dense_col);
      for (int d_row = 0; d_row < use_num_dense_col; d_row++) {
        int iCol = dense_col[d_row];
        double value = 0;
        for (int iEl = highs_a.start_[iCol]; iEl < highs_a.start_[iCol + 1];
             iEl++)
          value += hatA_d[offset + highs_a.index_[iEl]] * highs_a.value_[iEl];
        d_matrix[d_col][d_row] = value;
      }
      d_matrix[d_col][d_col] += 1 / theta_d[d_col];
      offset += system_size;
    }
    experiment_data.nla_time.factorization += getWallTime() - start_time;
  }
  invert.system_size = system_size;
  invert.use_num_dense_col = use_num_dense_col;
  experiment_data.nla_time.total = getWallTime() - start_time0;
  invert.valid = true;
  return kDecomposerStatusOk;
}

int newtonSolve(const HighsSparseMatrix &highs_a,
                const std::vector<double> &theta,
                const std::vector<double> &rhs, std::vector<double> &lhs,
                IpmInvert &invert, ExperimentData &experiment_data,
                const int& decomposer_source) {
  assert(invert.valid);

  double start_time = getWallTime();

  const int system_size = invert.system_size;
  const int use_num_dense_col = invert.use_num_dense_col;
  const std::vector<int> &dense_col = invert.dense_col;
  const std::vector<double> &theta_d = invert.theta_d;
  const std::vector<double> &hatA_d = invert.hatA_d;
  lhs = rhs;
  // Form \hat_b = \hat{A}_d^Tb (as d_rhs);
  if (decomposer_source == kDecomposerSourceSsids){
    callSsidsSolve(system_size, 1, lhs.data(), invert.ssids_data);
  } else if (decomposer_source == kDecomposerSourceMa86) {
    callMA86Solve(system_size, 1, lhs.data(), invert.ma86_data);
  } else if (decomposer_source == kDecomposerSourceQdldl){
    callQDLDLSolve(system_size, 1, lhs.data(), invert.qdldl_data);
  } else if (decomposer_source == kDecomposerSourceCholmod) {
    callCholmodSolve(system_size, 1, lhs.data(), invert.cholmod_data);
  } else if (decomposer_source == kDecomposerSourceHighs){
    callHighsSolve(system_size, 1, lhs.data(), invert.highs_factor);
  }
  if (use_num_dense_col) {
    std::vector<std::vector<double>> d_matrix = invert.d_matrix;
    std::vector<double> d_sol;
    std::vector<double> d_rhs;
    int offset = 0;
    for (int d_col = 0; d_col < use_num_dense_col; d_col++) {
      double value = 0;
      for (int iRow = 0; iRow < system_size; iRow++)
        value += hatA_d[offset + iRow] * rhs[iRow];
      d_rhs.push_back(value);
      offset += system_size;
    }
    // Now solve d_matrix.d_sol = d_rhs
    int gepp_status = gepp(d_matrix, d_rhs, d_sol);
    if (gepp_status)
      return kDecomposerStatusErrorFactorize;

    // Subtract \hat{A}_d^Td_sol from lhs
    offset = 0;
    for (int d_col = 0; d_col < use_num_dense_col; d_col++) {
      for (int iRow = 0; iRow < system_size; iRow++)
        lhs[iRow] -= hatA_d[offset + iRow] * d_sol[d_col];
      offset += system_size;
    }
  }
  experiment_data.nla_time.solve = getWallTime() - start_time;

  experiment_data.residual_error =
      residualErrorNewton(highs_a, theta, rhs, lhs);
  return kDecomposerStatusOk;
}

double newtonCondition(const HighsSparseMatrix &matrix,
                       const std::vector<double> &theta, IpmInvert &invert, const int& decomposer_source) {
  assert(invert.valid);
  SsidsData &ssids_data = invert.ssids_data;
  const bool chatty = false;

  // Get the largest eigenvalue by power method
  const double dim = matrix.num_row_;
  std::vector<double> from_iterate(dim);
  std::vector<double> to_iterate(dim);
  HighsRandom random;
  double lambda_n = 0;
  for (int iRow = 0; iRow < dim; iRow++)
    from_iterate[iRow] = random.fraction();
  for (int iter = 0; iter < 10; iter++) {
    double from_iterate_norm = Norm2(from_iterate);
    assert(from_iterate_norm > 0);
    VectorScale(from_iterate, 1 / from_iterate_norm);
    productAThetaAT(matrix, theta, from_iterate, to_iterate);
    double to_iterate_norm = Norm2(to_iterate);
    const bool converged =
        std::fabs(lambda_n - to_iterate_norm) / std::max(1.0, lambda_n) <
            1e-2 &&
        iter > 2;
    lambda_n = to_iterate_norm;
    if (chatty)
      printf("Iter %2d: lambda_n = %11.4g\n", iter, lambda_n);
    if (converged)
      break;
    from_iterate = to_iterate;
  }

  ExperimentData experiment_data;
  double lambda_1 = 0;
  for (int iRow = 0; iRow < dim; iRow++)
    from_iterate[iRow] = random.fraction();
  for (int iter = 0; iter < 10; iter++) {
    double from_iterate_norm = Norm2(from_iterate);
    assert(from_iterate_norm > 0);
    VectorScale(from_iterate, 1 / from_iterate_norm);
    newtonSolve(matrix, theta, from_iterate, to_iterate, invert,
                experiment_data, decomposer_source);
    double to_iterate_norm = Norm2(to_iterate);
    const bool converged =
        std::fabs(lambda_1 - to_iterate_norm) / std::max(1.0, lambda_1) <
            1e-2 &&
        iter > 2;
    lambda_1 = to_iterate_norm;
    if (chatty)
      printf("Iter %2d: lambda_1 = %11.4g\n", iter, lambda_1);
    if (converged)
      break;
    from_iterate = to_iterate;
  }
  double condition = lambda_n / lambda_1;
  if (chatty)
    printf("Condition = %g\n", condition);

  return condition;
}

bool increasingIndex(const HighsSparseMatrix &matrix) {
  if (matrix.isRowwise()) {
    for (int iRow = 0; iRow < matrix.num_row_; iRow++)
      for (int iEl = matrix.start_[iRow] + 1; iEl < matrix.start_[iRow + 1];
           iEl++)
        if (matrix.index_[iEl] <= matrix.index_[iEl - 1])
          return false;
  } else {
    for (int iCol = 0; iCol < matrix.num_col_; iCol++)
      for (int iEl = matrix.start_[iCol] + 1; iEl < matrix.start_[iCol + 1];
           iEl++)
        if (matrix.index_[iEl] <= matrix.index_[iEl - 1])
          return false;
  }
  return true;
}

void productAThetaAT(const HighsSparseMatrix &matrix,
                     const std::vector<double> &theta,
                     const std::vector<double> &x,
                     std::vector<double> &result) {
  assert(int(x.size()) == matrix.num_row_);
  std::vector<double> ATx;
  matrix.productTranspose(ATx, x);
  if (!theta.empty())
    for (int ix = 0; ix < matrix.num_col_; ix++)
      ATx[ix] *= theta[ix];
  matrix.product(result, ATx);
}

int computeAThetaAT(const HighsSparseMatrix &matrix,
		    const std::vector<double> &theta,
		    HighsSparseMatrix& AAT,
		    const int method,
		    const int max_num_nz) {
  // Create a row-wise copy of the matrix
  HighsSparseMatrix AT = matrix;
  AT.ensureRowwise();
  int num_zero_row = analyseScaledRowNorms(AT, theta, false);
  if (num_zero_row > 0) {
    printf("A.Theta has %d zero rows\n", num_zero_row);
    return kDecomposerStatusErrorFactorize;
  }

  int AAT_dim = matrix.num_row_;
  AAT.num_col_ = AAT_dim;
  AAT.num_row_ = AAT_dim;
  AAT.start_.resize(AAT_dim + 1, 0);

  std::vector<std::tuple<int, int, double>> non_zero_values;

  // First pass to calculate the number of non-zero elements in each column
  //
  if (method == 0) {
    std::vector<double> matrix_row(matrix.num_col_, 0);
    for (int iRow = 0; iRow < AAT_dim; iRow++) {
      for (int iEl = AT.start_[iRow]; iEl < AT.start_[iRow + 1]; iEl++)
        matrix_row[AT.index_[iEl]] = AT.value_[iEl];
      for (int iCol = iRow; iCol < AAT_dim; iCol++) {
        double dot = 0.0;
        for (int iEl = AT.start_[iCol]; iEl < AT.start_[iCol + 1]; iEl++) {
          const int ix = AT.index_[iEl];
          const double theta_i = !theta.empty() ? theta[ix] : 1;
          dot += theta_i * matrix_row[ix] * AT.value_[iEl];
        }
        if (dot != 0.0) {
          non_zero_values.emplace_back(iRow, iCol, dot);
          AAT.start_[iRow + 1]++;
          if (iRow != iCol)
            AAT.start_[iCol + 1]++;
        }
      }
      for (int iEl = AT.start_[iRow]; iEl < AT.start_[iRow + 1]; iEl++)
        matrix_row[AT.index_[iEl]] = 0;
      for (int ix = 0; ix < matrix.num_col_; ix++)
        assert(!matrix_row[ix]);
    }
  } else if (method == 1) {
    int AAT_num_nz = 0;
    std::vector<double> AAT_col_value(AAT_dim, 0);
    std::vector<int> AAT_col_index(AAT_dim);
    std::vector<bool> AAT_col_in_index(AAT_dim, false);
    for (int iRow = 0; iRow < AAT_dim; iRow++) {
      // Go along the row of A, and then down the columns corresponding
      // to its nonzeros
      int num_col_el = 0;
      for (int iRowEl = AT.start_[iRow]; iRowEl < AT.start_[iRow + 1]; iRowEl++) {
	int iCol = AT.index_[iRowEl];
	const double theta_value = !theta.empty() ? theta[iCol] : 1;
	if (!theta_value) continue;
	const double row_value = theta_value * AT.value_[iRowEl];
	for (int iColEl = matrix.start_[iCol]; iColEl < matrix.start_[iCol + 1]; iColEl++) {
	  int iRow1 = matrix.index_[iColEl];
	  if (iRow1 < iRow) continue;
	  double term = row_value * matrix.value_[iColEl];
	  if (!AAT_col_in_index[iRow1]) {
	    // This entry is not yet in the list of possible nonzeros
	    AAT_col_in_index[iRow1] = true;
	    AAT_col_index[num_col_el++] = iRow1;
	    AAT_col_value[iRow1] = term;
	  } else {
	    // This entry is in the list of possible nonzeros
	    AAT_col_value[iRow1] += term;
	  }
	}
      }
      for (int iEl = 0; iEl < num_col_el; iEl++) {
	int iCol = AAT_col_index[iEl];
	assert(iCol >= iRow);
	if (AAT_num_nz == max_num_nz) return kDecomposerStatusErrorOom;
	const double value = AAT_col_value[iCol];
	if (std::abs(value) > 1e-10) {
	  non_zero_values.emplace_back(iRow, iCol, value);
	  AAT.start_[iRow + 1]++;
	  if (iRow != iCol)
	    AAT.start_[iCol + 1]++;
	  AAT_num_nz += iRow != iCol ? 2 : 1;
	}
	AAT_col_value[iCol] = 0; // Not strictly necessary, but simplifies debugging
	AAT_col_in_index[iCol] = false;
      }
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
          AAT.start_[i + 1]++;
          if (i != j)
            AAT.start_[j + 1]++;
        }
      }
    }
  }

  // Prefix sum to get the correct column pointers
  for (int i = 0; i < AAT_dim; ++i)
    AAT.start_[i + 1] += AAT.start_[i];

  AAT.index_.resize(AAT.start_.back());
  AAT.value_.resize(AAT.start_.back());
  AAT.p_end_ = AAT.start_;
  AAT.p_end_.back() = AAT.index_.size();

  std::vector<int> current_positions = AAT.start_;

  // Second pass to actually fill in the indices and values
  for (const auto &val : non_zero_values) {
    int i = std::get<0>(val);
    int j = std::get<1>(val);
    double dot = std::get<2>(val);

    AAT.index_[current_positions[i]] = j;
    AAT.value_[current_positions[i]] = dot;
    current_positions[i]++;
    AAT.p_end_[i] = current_positions[i];

    if (i != j) {
      AAT.index_[current_positions[j]] = i;
      AAT.value_[current_positions[j]] = dot;
      current_positions[j]++;
      AAT.p_end_[j] = current_positions[j];
    }
  }
  AAT.p_end_.clear();
  return kDecomposerStatusOk;
}

// Gaussian elimination with partial pivoting for a dense matrix
//
// Returns 0 if it's successful
//
// Returns -1 if GE fails due to pivot less than kMinAbsPivot

const double kMinAbsPivot = 1e-12;
int gepp(const std::vector<std::vector<double>> &matrix,
         const std::vector<double> &rhs, std::vector<double> &solution) {
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
    if (pivot_row < 0)
      return -1;
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
    for (int iRow = k + 1; iRow < dim; iRow++) {
      double multiplier = ge_matrix[iRow][k] / ge_matrix[k][k];
      assert(std::fabs(multiplier) <= 1);
      for (int iCol = k + 1; iCol < dim; iCol++)
        ge_matrix[iRow][iCol] -= multiplier * ge_matrix[k][iCol];
      solution[iRow] -= multiplier * solution[k];
    }
  }
  for (int iRow = dim - 1; iRow >= 0; iRow--) {
    for (int iCol = iRow + 1; iCol < dim; iCol++)
      solution[iRow] -= solution[iCol] * ge_matrix[iRow][iCol];
    solution[iRow] /= ge_matrix[iRow][iRow];
  }
  return kDecomposerStatusOk;
}

int analyseScaledRowNorms(const HighsSparseMatrix &matrix,
			  const std::vector<double> &theta,
			  const bool quiet) {
  assert(matrix.isRowwise());
  int num_zero_row = 0;
  int num_row = matrix.num_row_;
  std::vector<double> row_norm(num_row, 0);
  for (int iRow = 0; iRow < num_row; iRow++) {
    double value = 0;
    for (int iEl = matrix.start_[iRow]; iEl < matrix.start_[iRow+1]; iEl++) 
      value += theta[matrix.index_[iEl]] * matrix.value_[iEl] * matrix.value_[iEl];
    row_norm[iRow] = value;
    if (value == 0) num_zero_row++;
  }
  if (!quiet) 
    analyseVectorValues(nullptr, "Squared 2-norm of A.Theta rows", matrix.num_row_,
			row_norm);
  return num_zero_row;
}

// int ssids_decompose();

int callSsidsAugmentedFactor(const HighsSparseMatrix &matrix,
                             const std::vector<double> &theta,
                             SsidsData &ssids_data,
                             ExperimentData &experiment_data) {
#ifdef HAVE_SPRAL
  experiment_data.nla_time.form = 0;
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
    val.push_back(-1 / theta_i);
    row.push_back(iCol + array_base);
    for (int iEl = matrix.start_[iCol]; iEl < matrix.start_[iCol + 1]; iEl++) {
      int iRow = matrix.index_[iEl];
      val.push_back(matrix.value_[iEl]);
      row.push_back(row_index_offset + iRow + array_base);
    }
  }
  // Now the (zero) columns for [A^T]
  //                            [ 0 ]

  const double diagonal = ssids_data.options.small; // 1e-20
  for (int iRow = 0; iRow < matrix.num_row_; iRow++) {
    ptr.push_back(val.size());
    if (diagonal) {
      val.push_back(diagonal);
      row.push_back(row_index_offset + iRow + array_base);
    }
  }
  // Possibly add 1 to all the pointers
  if (array_base)
    for (auto &p : ptr)
      ++p;

  // Add the last pointer
  ptr.push_back(val.size() + array_base);

  long *ptr_ptr = ptr.data();
  int *row_ptr = row.data();
  double *val_ptr = val.data();

  // Initialize derived types
  ssids_data.akeep = nullptr;
  ssids_data.fkeep = nullptr;
  // Need to set to 1 if using Fortran 1-based indexing
  ssids_data.options.array_base = array_base;
  //  ssids_data.options.print_level = 2;

  experiment_data.system_nnz = matrix.num_col_ + 2 * matrix.numNz() + matrix.num_row_;
  experiment_data.nla_time.setup = getWallTime() - start_time;

  // Perform analyse and factorise with data checking
  bool check = true;
  start_time = getWallTime();
  spral_ssids_analyse(check, experiment_data.system_size, nullptr, ptr_ptr,
		      row_ptr, nullptr, &ssids_data.akeep, &ssids_data.options,
		      &ssids_data.inform);
  experiment_data.nla_time.analysis = getWallTime() - start_time;
  if (ssids_data.inform.flag < 0) {
    printf("SSIDS analyse:   flag = %d\n", ssids_data.inform.flag);
    return kDecomposerStatusErrorFactorize;
  }
  bool positive_definite = false;
  start_time = getWallTime();
  spral_ssids_factor(positive_definite, nullptr, nullptr, val_ptr, nullptr,
		     ssids_data.akeep, &ssids_data.fkeep, &ssids_data.options,
		     &ssids_data.inform);
  experiment_data.nla_time.factorization = getWallTime() - start_time;
  if (ssids_data.inform.flag < 0) {
    printf("SSIDS factorize: flag = %d\n", ssids_data.inform.flag);
    return kDecomposerStatusErrorFactorize;
  }
  experiment_data.nnz_L = ssids_data.inform.num_factor;
  experiment_data.fillIn_LL();
  return kDecomposerStatusOk;
#else
  return kDecomposerStatusErrorFactorize;
#endif
}

int callSsidsNewtonFactor(const HighsSparseMatrix &AThetaAT,
                          SsidsData &ssids_data,
                          ExperimentData &experiment_data) {
#ifdef HAVE_SPRAL
  double start_time = getWallTime();
  // Prepare data structures for SPRAL
  const int array_base = 0;
  std::vector<long> ptr;
  std::vector<int> row;
  std::vector<double> val;
  spral_ssids_default_options(&ssids_data.options);

  // Extract lower triangular part of AAT
  for (int col = 0; col < AThetaAT.num_col_; col++) {
    ptr.push_back(val.size());
    for (int idx = AThetaAT.start_[col]; idx < AThetaAT.start_[col + 1];
         idx++) {
      int row_idx = AThetaAT.index_[idx];
      if (row_idx >= col) {
        val.push_back(AThetaAT.value_[idx]);
        row.push_back(row_idx + array_base);
      }
    }
  }
  // Possibly add 1 to all the pointers
  if (array_base)
    for (auto &p : ptr)
      ++p;

  // Add the last pointer
  ptr.push_back(val.size() + array_base);

  long *ptr_ptr = ptr.data();
  int *row_ptr = row.data();
  double *val_ptr = val.data();

  // Initialize derived types
  ssids_data.akeep = nullptr;
  ssids_data.fkeep = nullptr;
  // Need to set to 1 if using Fortran 1-based indexing
  ssids_data.options.array_base =
    array_base; // Need to set to 1 if using Fortran 1-based indexing

  experiment_data.nla_time.setup = getWallTime() - start_time;

  // Perform analyse and factorise with data checking
  bool check = true;
  start_time = getWallTime();
  spral_ssids_analyse(check, AThetaAT.num_col_, nullptr, ptr_ptr, row_ptr,
		      nullptr, &ssids_data.akeep, &ssids_data.options,
		      &ssids_data.inform);
  experiment_data.nla_time.analysis = getWallTime() - start_time;
  if (ssids_data.inform.flag < 0) {
    printf("SSIDS analyse:   flag = %d\n", ssids_data.inform.flag);
    return kDecomposerStatusErrorFactorize;
  }

  bool positive_definite = true;
  start_time = getWallTime();
  spral_ssids_factor(positive_definite, nullptr, nullptr, val_ptr, nullptr,
		     ssids_data.akeep, &ssids_data.fkeep, &ssids_data.options,
		     &ssids_data.inform);
  experiment_data.nla_time.factorization = getWallTime() - start_time;
  if (ssids_data.inform.flag < 0) {
    printf("SSIDS factorize: flag = %d\n", ssids_data.inform.flag);
    return kDecomposerStatusErrorFactorize;
  }
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
  return kDecomposerStatusOk;
#else
  return kDecomposerStatusErrorFactorize;
#endif
}

void callSsidsSolve(const int system_size, const int num_rhs, double *rhs,
                    SsidsData &ssids_data) {
#ifdef HAVE_SPRAL
  spral_ssids_solve(0, num_rhs, rhs, system_size, ssids_data.akeep,
		    ssids_data.fkeep, &ssids_data.options, &ssids_data.inform);
#endif
}

int callMA86AugmentedFactor(const HighsSparseMatrix& matrix,
			    const std::vector<double>& theta,
			    MA86Data& ma86_data,
			    ExperimentData& experiment_data){
#ifdef HAVE_MA86
  int n = matrix.num_col_;
  int m = matrix.num_row_;
  experiment_data.nla_time.form = 0;
  double start_time = getWallTime();

  const int array_base = 0;
  std::vector<int> ptr;
  std::vector<int> row;
  std::vector<double> val;

  // Extract lower triangular part of augmented system
  int row_index_offset = matrix.num_col_;
  for (int iCol = 0; iCol < matrix.num_col_; iCol++) {
    const double theta_i = !theta.empty() ? theta[iCol] : 1;
    
    ptr.push_back(val.size()); 
    val.push_back(-1/theta_i);
    row.push_back(iCol); 
    for (int iEl = matrix.start_[iCol]; iEl < matrix.start_[iCol+1]; iEl++) {
      int iRow = matrix.index_[iEl];

      val.push_back(matrix.value_[iEl]);
      row.push_back(row_index_offset + iRow);
    }
  }

  const double diagonal = 1e-20;//1e-20;
  for (int iRow = 0; iRow < matrix.num_row_; iRow++) {
    ptr.push_back(val.size()); 
    if (diagonal) {
      val.push_back(diagonal);
      row.push_back(row_index_offset + iRow); 
    }
  }

  ptr.push_back(val.size());

  ma86_data.order.resize(m+n);
  for (int i = 0; i < m+n; i++) ma86_data.order[i] = i;

  wrapper_ma86_default_control(&ma86_data.control);

  experiment_data.system_nnz = matrix.num_col_ + 2 * highs_a.numNz();// + matrix.num_row;
  experiment_data.setup_time = getWallTime() - start_time;

  start_time = getWallTime();
  wrapper_ma86_analyse(n+m, ptr.data(), row.data(), ma86_data.order.data(),
                       &ma86_data.keep, &ma86_data.control, &ma86_data.info);
  experiment_data.analysis_time = getWallTime() - start_time;
  if (ma86_data.info.flag < 0) return kDecomposerStatusErrorFactorize;
  start_time = getWallTime();
  wrapper_ma86_factor(n+m, ptr.data(), row.data(), val.data(), ma86_data.order.data(),
                       &ma86_data.keep, &ma86_data.control, &ma86_data.info);
  if (ma86_data.info.flag < 0) return kDecomposerStatusErrorFactorize;
  experiment_data.factorization_time = getWallTime() - start_time;
  experiment_data.nnz_L = ma86_data.info.num_factor;
  experiment_data.fillIn_LL();

  return kDecomposerStatusOk;
#else
  return kDecomposerStatusErrorFactorize;
#endif
}

int callMA86NewtonFactor(const HighsSparseMatrix& AThetaAT,
			                  MA86Data& ma86_data,
			                  ExperimentData& experiment_data){
#ifdef HAVE_MA86
  double start_time = getWallTime();

  int n = AThetaAT.num_col_;
  const int array_base = 0;
  std::vector<int> ptr;
  std::vector<int> row;
  std::vector<double> val;

  // Extract upper triangular part of AAT in fortran
  for (int col = 0; col < AThetaAT.num_col_; col++){
      ptr.push_back(val.size());
      for (int idx = AThetaAT.start_[col]; idx < AThetaAT.start_[col+1]; idx++){
          int row_idx = AThetaAT.index_[idx];
          if (row_idx >= col){
              val.push_back(AThetaAT.value_[idx]);
              row.push_back(row_idx);
          }
      }
  }

  // Add the last pointer
  ptr.push_back(val.size());

  ma86_data.order.resize(AThetaAT.num_row_);
  for (int i = 0; i < n; i++) ma86_data.order[i] = i;

  wrapper_ma86_default_control(&ma86_data.control);
  experiment_data.setup_time = getWallTime() - start_time;
  start_time = getWallTime();
  wrapper_ma86_analyse(n, ptr.data(), row.data(), ma86_data.order.data(),
                       &ma86_data.keep, &ma86_data.control, &ma86_data.info);
  experiment_data.analysis_time = getWallTime() - start_time;
  if (ma86_data.info.flag < 0) return kDecomposerStatusErrorFactorize;
  start_time = getWallTime();
  wrapper_ma86_factor(n, ptr.data(), row.data(), val.data(), ma86_data.order.data(),
                       &ma86_data.keep, &ma86_data.control, &ma86_data.info);
  if (ma86_data.info.flag < 0) return kDecomposerStatusErrorFactorize;
  experiment_data.factorization_time = getWallTime() - start_time;
  experiment_data.nnz_L = ma86_data.info.num_factor;
  experiment_data.fillIn_LL();

  return kDecomposerStatusOk;
#else
  return kDecomposerStatusErrorFactorize;
#endif
}

void callMA86Solve(const int system_size,
		   const int num_rhs,
		   double* rhs,
		   MA86Data& ma86_data){
#ifdef HAVE_MA86
  wrapper_ma86_solve(0,1, system_size, rhs, ma86_data.order.data(), &ma86_data.keep, &ma86_data.control, &ma86_data.info);
#endif
}

int callQDLDLNewtonFactor(const HighsSparseMatrix& AThetaAT,
			 QDLDLData& qdldl_data,
			 ExperimentData& experiment_data){
#ifdef HAVE_QDLDL
  double start_time = getWallTime();

  int n = AThetaAT.num_col_;
  QDLDL_int An = AThetaAT.num_col_;
  std::vector<QDLDL_int> Ap(An+1);
  std::vector<QDLDL_int> Ai;
  std::vector<QDLDL_float> Ax;

  for (int col = 0; col < AThetaAT.num_col_; col++){
      Ap[col] = Ai.size();
      for (int idx = AThetaAT.start_[col]; idx < AThetaAT.start_[col+1]; idx++){
          int row_idx = AThetaAT.index_[idx];
          if (row_idx <= col){
              Ai.push_back(row_idx);
              Ax.push_back(AThetaAT.value_[idx]);
          }
      }
  }
  Ap[An] = Ai.size();

  QDLDL_int *Ap_star = Ap.data();
  QDLDL_int *Ai_star = Ai.data();
  QDLDL_float *Ax_star = Ax.data();
 
  experiment_data.setup_time = getWallTime() - start_time;
  /*--------------------------------
  * pre-factorisation memory allocations
  *---------------------------------*/
  start_time = getWallTime();
  qdldl_data.etree = (QDLDL_int*)malloc(sizeof(QDLDL_int)*An);
  qdldl_data.Lnz   = (QDLDL_int*)malloc(sizeof(QDLDL_int)*An);

  //For the L factors.   Li and Lx are sparsity dependent
  //so must be done after the etree is constructed
  qdldl_data.Lp    = (QDLDL_int*)malloc(sizeof(QDLDL_int)*(An+1));
  qdldl_data.D     = (QDLDL_float*)malloc(sizeof(QDLDL_float)*An);
  qdldl_data.Dinv  = (QDLDL_float*)malloc(sizeof(QDLDL_float)*An);


  qdldl_data.iwork = (QDLDL_int*)malloc(sizeof(QDLDL_int)*(3*An));
  qdldl_data.bwork = (QDLDL_bool*)malloc(sizeof(QDLDL_bool)*An);
  qdldl_data.fwork = (QDLDL_float*)malloc(sizeof(QDLDL_float)*An);

  
  /*--------------------------------
  * elimination tree calculation
  *---------------------------------*/
  
  qdldl_data.sumLnz = QDLDL_etree(An,Ap_star,Ai.data(),qdldl_data.iwork,qdldl_data.Lnz,qdldl_data.etree);
  experiment_data.analysis_time = getWallTime() - start_time;
  /*--------------------------------
  * LDL factorisation
  *---------------------------------*/
  start_time = getWallTime();
  //First allocate memory for Li and Lx
  experiment_data.nnz_L = qdldl_data.sumLnz;
  //std::cout << "sumLnz = " << sumLnz << std::endl;
  qdldl_data.Li    = (QDLDL_int*)malloc(sizeof(QDLDL_int)*qdldl_data.sumLnz);
  
  qdldl_data.Lx    = (QDLDL_float*)malloc(sizeof(QDLDL_float)*qdldl_data.sumLnz);

  //now factor
  QDLDL_factor(An,Ap_star,Ai.data(),Ax.data(),qdldl_data.Lp,qdldl_data.Li,qdldl_data.Lx,qdldl_data.D,qdldl_data.Dinv,qdldl_data.Lnz,qdldl_data.etree,qdldl_data.bwork,qdldl_data.iwork,qdldl_data.fwork);
  qdldl_data.x = (QDLDL_float*)malloc(sizeof(QDLDL_float)*An);
  experiment_data.factorization_time = getWallTime() - start_time;
  experiment_data.fillIn_LDL();
  return kDecomposerStatusOk;
#else
  return kDecomposerStatusErrorFactorize;
#endif
}
/*
void callQDLDLSolve(const int system_size,
		   const int num_rhs,
		   double* rhs,
		   QDLDLData& qdldl_data){
  // print rhs
  QDLDL_solve(qdldl_data.Ln,qdldl_data.Lp,qdldl_data.Li,qdldl_data.Lx,qdldl_data.Dinv, qdldl_data.x);
      
}*/
void callQDLDLSolve(const int system_size,
                    const int num_rhs,
                    double* rhs,
                    QDLDLData& qdldl_data){
#ifdef HAVE_QDLDL
  // Copy rhs to qdldl_data.x before the solve operation
  //memcpy(qdldl_data.x, rhs, system_size * sizeof(double));

  //for(int i=0;i < system_size-1; i++) qdldl_data.x[i] = rhs[i];
  // Solve the system
  QDLDL_solve(qdldl_data.Ln,qdldl_data.Lp,qdldl_data.Li,qdldl_data.Lx,qdldl_data.Dinv, rhs);
  
  // Copy the result back into rhs
  //memcpy(rhs, qdldl_data.x, system_size * sizeof(double));
#endif
}


int callQDLDLAugmentedFactor(const HighsSparseMatrix& matrix,
			    const std::vector<double>& theta,
			    QDLDLData& qdldl_data,
			    ExperimentData& experiment_data){
#ifdef HAVE_QDLDL
  experiment_data.nla_time.form = 0;          
  double start_time = getWallTime();
  HighsSparseMatrix A_rowwise=matrix;
  A_rowwise.ensureRowwise();

  int n = matrix.num_col_;
  QDLDL_int An = matrix.num_col_;
  const int array_base = 0;
  std::vector<QDLDL_int> Ap;
  std::vector<QDLDL_int> Ai;
  std::vector<QDLDL_float> Ax;

  int row_index_offset = matrix.num_col_;
  int column_index_offset = matrix.num_row_;
  // Extract lower triangular part of augmented system row-wise
  //
  // First the rows    for [-\Theta^{-1}]
  //                       [      A     ]
  for (int iCol = 0; iCol < A_rowwise.num_col_; iCol++){
    const double theta_i = !theta.empty() ? theta[iCol] : 1;
    Ap.push_back(Ax.size());
    Ax.push_back(-1/theta_i);
    Ai.push_back(iCol + array_base);
    
  }

  const double diagonal = 0; //1e-20
  for (int iRow = 0; iRow < matrix.num_row_; iRow++){
    Ap.push_back(Ax.size());
    for (int idx = A_rowwise.start_[iRow]; idx < A_rowwise.start_[iRow+1]; idx++){
      int row_idx = A_rowwise.index_[idx];
      Ai.push_back(row_idx + array_base );
      Ax.push_back(A_rowwise.value_[idx]);  
    }
    
    if (diagonal) {
      Ai.push_back(iRow + row_index_offset + array_base);
      Ax.push_back(diagonal);
    }
  }

  Ap.push_back(Ax.size());

  QDLDL_int *Ap_star = Ap.data();
  QDLDL_int *Ai_star = Ai.data();
  QDLDL_float *Ax_star = Ax.data();
  experiment_data.system_nnz = matrix.num_col_ + 2 * highs_a.numNz()
  experiment_data.setup_time = getWallTime() - start_time;
  /*--------------------------------
  * pre-factorisation memory allocations
  *---------------------------------*/
  start_time = getWallTime();
  qdldl_data.etree = (QDLDL_int*)malloc(sizeof(QDLDL_int)*An);
  qdldl_data.Lnz   = (QDLDL_int*)malloc(sizeof(QDLDL_int)*An);

  //For the L factors.   Li and Lx are sparsity dependent
  //so must be done after the etree is constructed
  qdldl_data.Lp    = (QDLDL_int*)malloc(sizeof(QDLDL_int)*(An+1));
  qdldl_data.D     = (QDLDL_float*)malloc(sizeof(QDLDL_float)*An);
  qdldl_data.Dinv  = (QDLDL_float*)malloc(sizeof(QDLDL_float)*An);


  qdldl_data.iwork = (QDLDL_int*)malloc(sizeof(QDLDL_int)*(3*An));
  qdldl_data.bwork = (QDLDL_bool*)malloc(sizeof(QDLDL_bool)*An);
  qdldl_data.fwork = (QDLDL_float*)malloc(sizeof(QDLDL_float)*An);

  
  /*--------------------------------
  * elimination tree calculation
  *---------------------------------*/
  
  qdldl_data.sumLnz = QDLDL_etree(An,Ap_star,Ai.data(),qdldl_data.iwork,qdldl_data.Lnz,qdldl_data.etree);
  experiment_data.analysis_time = getWallTime() - start_time;
  /*--------------------------------
  * LDL factorisation
  *---------------------------------*/
  start_time = getWallTime();
  //First allocate memory for Li and Lx
  experiment_data.nnz_L = qdldl_data.sumLnz;
  std::cout << "sumLnz = " << qdldl_data.sumLnz << std::endl;
  qdldl_data.Li    = (QDLDL_int*)malloc(sizeof(QDLDL_int)*qdldl_data.sumLnz);
  
  qdldl_data.Lx    = (QDLDL_float*)malloc(sizeof(QDLDL_float)*qdldl_data.sumLnz);
  qdldl_data.x = (QDLDL_float*)malloc(sizeof(QDLDL_float)*An);
  
  //now factor
  QDLDL_factor(An,Ap_star,Ai.data(),Ax.data(),qdldl_data.Lp,qdldl_data.Li,qdldl_data.Lx,qdldl_data.D,qdldl_data.Dinv,qdldl_data.Lnz,qdldl_data.etree,qdldl_data.bwork,qdldl_data.iwork,qdldl_data.fwork);
  experiment_data.factorization_time = getWallTime() - start_time;
  experiment_data.fillIn_LDL();
  return kDecomposerStatusOk;
#else
  return kDecomposerStatusErrorFactorize;
#endif
}

int callCholmodNewtonFactor(const HighsSparseMatrix &AThetaAT,
                          CholmodData &cholmod_data,
                          ExperimentData &experiment_data,
                          const int num_thods){
#ifdef HAVE_CHOLMOD
  experiment_data.nla_time.form = 0;
  double start_time = getWallTime();

  cholmod_data.c;
  cholmod_start(&cholmod_data.c);

  //create the cholmod_triplet structure
  //cholmod_triplet* T = cholmod_allocate_triplet(AThetaAT.num_row_, AThetaAT.num_col_, AThetaAT.value_.size(), 1, CHOLMOD_REAL, &cholmod_data.c);
  cholmod_data.T = cholmod_allocate_triplet(AThetaAT.num_row_, AThetaAT.num_col_, AThetaAT.value_.size(), 1, CHOLMOD_REAL, &cholmod_data.c);
  int* Ti = static_cast<int*>(cholmod_data.T->i);
  int* Tj = static_cast<int*>(cholmod_data.T->j);
  double* Tx = static_cast<double*>(cholmod_data.T->x);
  // Extract upper triangular part of AThetaAT
  for (int col = 0; col < AThetaAT.num_col_; col++){
      for (int idx = AThetaAT.start_[col]; idx < AThetaAT.start_[col+1]; idx++){
          if (col >= AThetaAT.index_[idx]){
              *Ti++ = AThetaAT.index_[idx];
              *Tj++ = col;
              *Tx++ = AThetaAT.value_[idx];
          }
      }
  }
  cholmod_data.T->nnz = Tx - static_cast<double*>(cholmod_data.T->x);

//   cholmod_print_triplet(T, "T", &c);
  cholmod_check_triplet(cholmod_data.T, &cholmod_data.c);
  
  cholmod_data.a = cholmod_triplet_to_sparse(cholmod_data.T, 1, &cholmod_data.c);

//    cholmod_print_sparse(a,"A",&c);
  cholmod_check_sparse(cholmod_data.a,&cholmod_data.c);
  experiment_data.setup_time = getWallTime() - start_time;
  start_time = getWallTime();
  /*
  If you are going to factorize hundreds or more matrices with the same
  * nonzero pattern, you may wish to spend a great deal of time finding a
  * good permutation.  In this case, try setting Common->nmethods to 9.
  * The time spent in cholmod_analysis will be very high, but you need to
  * call it only once.
  */
  //c.nmethods = num_thods;
  //c.method[2].ordering = CHOLMOD_METIS;
  //c.postorder = 0;
  cholmod_data.L = cholmod_analyze(cholmod_data.a, &cholmod_data.c);
  experiment_data.analysis_time = getWallTime() - start_time;
  start_time = getWallTime();

  cholmod_factorize(cholmod_data.a, cholmod_data.L, &cholmod_data.c);
  
  size_t nnz = 0;
  if (cholmod_data.L->is_super) {  // Supernodal factor
      for (size_t s = 0; s < cholmod_data.L->nsuper; s++) {
          int pi = ((int*)(cholmod_data.L->pi))[s];
          int pj = ((int*)(cholmod_data.L->pi))[s+1];
          for (int p = pi; p < pj; p++) {
              int i = ((int*)(cholmod_data.L->s))[p];
              nnz += (cholmod_data.L->n - i);
          }
      }
      //std::cout << "Number of nonzeros in L = " << nnz << std::endl;
  }
  else {  // Simplicial factor
      for (size_t j = 0; j < cholmod_data.L->n; j++) {
          nnz += ((int*)(cholmod_data.L->nz))[j];
      }
      //std::cout << "Number of nonzeros in L = " << nnz << std::endl;
  }
  experiment_data.factorization_time = getWallTime() - start_time;
  experiment_data.nnz_L = nnz;
  experiment_data.fillIn_LL();

  return kDecomposerStatusOk;
#else
  return kDecomposerStatusErrorFactorize;
#endif
}

void callCholmodSolve(const int system_size, const int num_rhs, double *rhs,
                    CholmodData &cholmod_data){
#ifdef HAVE_CHOLMOD
  // Create a cholmod_dense structure for rhs
  cholmod_dense* rhs_dense = cholmod_allocate_dense(system_size, num_rhs, system_size, CHOLMOD_REAL, &(cholmod_data.c));
  std::copy_n(rhs, system_size, static_cast<double*>(rhs_dense->x));

  cholmod_dense* x = cholmod_solve(CHOLMOD_A, cholmod_data.L, rhs_dense, &(cholmod_data.c));

  // Copy the solution back into rhs
  std::copy_n(static_cast<double*>(x->x), system_size, rhs);

  cholmod_free_dense(&rhs_dense, &(cholmod_data.c));
#endif
}

int callHighsAugmentedFactor(const HighsSparseMatrix &matrix,
                             const std::vector<double> &theta,
                             HFactor& highs_factor,
                             ExperimentData &experiment_data) {
  experiment_data.nla_time.form = 0;
  double start_time = getWallTime();

  // Prepare data structures for HiGHS
  HighsSparseMatrix augmented;
  // Extract lower triangular part of augmented system
  //
  // First the columns for [-\Theta^{-1}]
  //                       [      A     ]
  int row_index_offset = matrix.num_col_;
  std::vector<HighsInt>& start = augmented.start_;
  std::vector<HighsInt>& index = augmented.index_;
  std::vector<double>& value = augmented.value_;
  for (int iCol = 0; iCol < matrix.num_col_; iCol++) {
    const double theta_i = !theta.empty() ? theta[iCol] : 1;
    value.push_back(-1 / theta_i);
    index.push_back(iCol);
    for (int iEl = matrix.start_[iCol]; iEl < matrix.start_[iCol + 1]; iEl++) {
      value.push_back(matrix.value_[iEl]);
      index.push_back(row_index_offset + matrix.index_[iEl]);
    }
    start.push_back(index.size());
  }
  // Now the columns for [A^T]
  //                     [ 0 ]
  for (int iRow = 0; iRow < matrix.num_row_; iRow++) {
  // !!!! Insert A^T efficiently
    start.push_back(index.size());
  }
  
  return kDecomposerStatusErrorFactorize;
}

int callHighsNewtonFactor(const HighsSparseMatrix &AThetaAT,
			  HFactor& highs_factor,
			  ExperimentData &experiment_data) {
  return kDecomposerStatusErrorFactorize;
}

void callHighsSolve(const int system_size, const int num_rhs, double *rhs,
		    HFactor& highs_factor) {
}

