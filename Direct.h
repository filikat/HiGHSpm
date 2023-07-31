#ifndef DIRECT_H
#define DIRECT_H

#include "spral.h"
#include "ExperimentData.h"
#include "util/HighsSparseMatrix.h"

struct SsidsData {
  void *akeep{nullptr};
  void *fkeep{nullptr};
  struct spral_ssids_options options;
  struct spral_ssids_inform inform;
  int clear();
};

struct IpmInvert {
  bool valid = false;
  int system_size = 0;
  int use_num_dense_col = 0;
  std::vector<int> dense_col;
  std::vector<double> theta_d;
  std::vector<double> hatA;
  SsidsData ssids_data;
  int clear();
};

int augmentedInvert(const HighsSparseMatrix &highs_a,
		   const std::vector<double> &theta,
		   IpmInvert& invert,
		   ExperimentData& experiment_data);

void augmentedSolve(const HighsSparseMatrix &highs_a,
		   const std::vector<double> &theta,
		   const std::vector<double> &rhs_x,
		   const std::vector<double> &rhs_y,
		   std::vector<double> &lhs_x,
		   std::vector<double> &lhs_y,
		   IpmInvert& invert,
		   ExperimentData& experiment_data);

int newtonInvert(const HighsSparseMatrix &highs_a,
		 const std::vector<double> &theta,
		 IpmInvert& invert,
		 const int option_max_dense_col,
		 const double option_dense_col_tolerance,
		 ExperimentData& experiment_data);

int newtonSystemSolve(const HighsSparseMatrix &highs_a,
		const std::vector<double> &theta,
		const std::vector<double> &rhs,
		std::vector<double> &lhs,
		IpmInvert& invert,
		ExperimentData& experiment_data);

int newtonSolve(const HighsSparseMatrix &highs_a,
		const std::vector<double> &theta,
		const std::vector<double> &rhs,
		std::vector<double> &lhs,
		IpmInvert& invert,
		const int option_max_dense_col,
		const double option_dense_col_tolerance,
		ExperimentData& experiment_data);

bool increasingIndex(const HighsSparseMatrix& matrix);
void productAThetaAT(const HighsSparseMatrix& matrix,
		     const std::vector<double>& theta,
		     const std::vector<double>& x,
		     std::vector<double>& result);
HighsSparseMatrix computeAThetaAT(const HighsSparseMatrix& matrix,
				  const std::vector<double>& theta);

int gepp(const std::vector<std::vector<double>>& matrix,
	 const std::vector<double>& rhs,
	 std::vector<double>& solution);

int callSsidsAugmentedFactor(const HighsSparseMatrix& matrix,
			     const std::vector<double>& theta,
			     SsidsData& ssids_data,
			     ExperimentData& experiment_data);
int callSsidsNewtonFactor(const HighsSparseMatrix& AThetaAT,
			  SsidsData& ssids_data,
			  ExperimentData& experiment_data);
void callSsidsSolve(const int system_size,
		    const int num_rhs,
		    double* rhs,
		    SsidsData& ssids_data);
#endif
