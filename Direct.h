#ifndef DIRECT_H
#define DIRECT_H

#include "ExperimentData.h"
#include "VectorOperations.h"
#include "spral.h" 
#include "util/HighsSparseMatrix.h"
#include "suitesparse/cholmod.h"
extern "C" {
	#include "hsl_ma86_wrapper.h"
  #include "qdldl.h"
}

struct SsidsData {
  void *akeep{nullptr};
  void *fkeep{nullptr};
  struct spral_ssids_options options;
  struct spral_ssids_inform inform;
  int clear();
};

struct MA86Data { 
	void *keep;
	ma86_control_d control;
	ma86_info_d info;
	std::vector<int> order;
	void clear();
};

struct QDLDLData{
  QDLDL_int Ln;
  QDLDL_int *Lp;
  QDLDL_int *Li;
  QDLDL_float *Lx;
  QDLDL_float *D;
  QDLDL_float *Dinv;

  QDLDL_int *etree;
  QDLDL_int *Lnz;
  QDLDL_int sumLnz;

  QDLDL_int   *iwork;
  QDLDL_bool  *bwork;
  QDLDL_float *fwork;

  QDLDL_float *x; //Data for results of A\b
  void clear();
};

struct CholmodData{
  cholmod_common c;
  cholmod_triplet* T;
  cholmod_sparse* a;
  cholmod_factor* L;
  cholmod_dense* x;
  cholmod_dense* b;
  void clear();
};

struct IpmInvert {
  bool valid = false;
  int system_size = 0;
  int use_num_dense_col = 0;
  std::vector<int> dense_col;
  std::vector<double> theta_d;
  std::vector<double> hatA_d;
  std::vector<std::vector<double>> d_matrix;
  SsidsData ssids_data;
  MA86Data ma86_data;
  QDLDLData qdldl_data;
  CholmodData cholmod_data;
  void clear(const int solver_type = 1);
};

void chooseDenseColumns(const HighsSparseMatrix &highs_a,
			const std::vector<double> &theta, 
			const int option_max_dense_col,
			const double option_dense_col_tolerance,
			std::vector<int> &dense_col,
			ExperimentData &experiment_data,
			const bool quiet = true);

int augmentedInvert(const HighsSparseMatrix &highs_a,
                    const std::vector<double> &theta, IpmInvert &invert,
                    ExperimentData &experiment_data,
		                const int solver_type = 1);

void augmentedSolve(const HighsSparseMatrix &highs_a,
                    const std::vector<double> &theta,
                    const std::vector<double> &rhs_x,
                    const std::vector<double> &rhs_y,
                    std::vector<double> &lhs_x, std::vector<double> &lhs_y,
                    IpmInvert &invert, ExperimentData &experiment_data,
                    const int solver_type = 1);

double augmentedCondition(const HighsSparseMatrix &matrix,
                          const std::vector<double> &theta, IpmInvert &invert,
                          const int solver_type = 1);

int newtonInvert(const HighsSparseMatrix &highs_a,
		 const std::vector<double> &theta,
		 IpmInvert& invert,
		 const int option_max_dense_col,
		 const double option_dense_col_tolerance,
		 ExperimentData& experiment_data,
		 const bool quiet = true,
		 const int solver_type = 1);

int newtonSolve(const HighsSparseMatrix &highs_a,
                const std::vector<double> &theta,
                const std::vector<double> &rhs, std::vector<double> &lhs,
                IpmInvert &invert, ExperimentData &experiment_data,
                const int &solver_type = 1);

double newtonCondition(const HighsSparseMatrix &matrix,
                       const std::vector<double> &theta, IpmInvert &invert, 
                       const int &solver_type = 1);

bool increasingIndex(const HighsSparseMatrix &matrix);
void productAThetaAT(const HighsSparseMatrix &matrix,
                     const std::vector<double> &theta,
                     const std::vector<double> &x, std::vector<double> &result);
HighsSparseMatrix computeAThetaAT(const HighsSparseMatrix &matrix,
                                  const std::vector<double> &theta);

int gepp(const std::vector<std::vector<double>> &matrix,
         const std::vector<double> &rhs, std::vector<double> &solution);

int callSsidsAugmentedFactor(const HighsSparseMatrix &matrix,
                             const std::vector<double> &theta,
                             SsidsData &ssids_data,
                             ExperimentData &experiment_data);
int callSsidsNewtonFactor(const HighsSparseMatrix &AThetaAT,
                          SsidsData &ssids_data,
                          ExperimentData &experiment_data);
void callSsidsSolve(const int system_size, const int num_rhs, double *rhs,
                    SsidsData &ssids_data);

int callMA86AugmentedFactor(const HighsSparseMatrix& matrix,
			    const std::vector<double>& theta,
			    MA86Data& ma86_data,
			    ExperimentData& experiment_data);

int callMA86NewtonFactor(const HighsSparseMatrix& AThetaAT,
			 MA86Data& ma86_data,
			 ExperimentData& experiment_data);

void callMA86Solve(const int system_size,
		   const int num_rhs,
		   double* rhs,
		   MA86Data& ma86_data);		

int callQDLDLAugmentedFactor(const HighsSparseMatrix& matrix,
			    const std::vector<double>& theta,
			    QDLDLData& qdldl_data,
			    ExperimentData& experiment_data);

int callQDLDLNewtonFactor(const HighsSparseMatrix& AThetaAT,
			 QDLDLData& qdldl_data,
			 ExperimentData& experiment_data);

void callQDLDLSolve(const int system_size,
		   const int num_rhs,
		   double* rhs,
		   QDLDLData& qdldl_data);	

int callCholmodNewtonFactor(const HighsSparseMatrix &AThetaAT,
                          CholmodData &cholmod_data,
                          ExperimentData &experiment_data,
                          const int num_thods = 3);
void callCholmodSolve(const int system_size, const int num_rhs, double *rhs,
                    CholmodData &cholmod_data);	

#endif
