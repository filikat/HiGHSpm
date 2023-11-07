#ifndef DIRECT_H
#define DIRECT_H

#include "ExperimentData.h"
#include "IPM_const.h"
#include "VectorOperations.h"
#ifdef HAVE_SPRAL
#include "spral.h"
#endif
#include "util/HighsSparseMatrix.h"
#ifdef HAVE_CHOLMOD
#include "cholmod.h"
#endif
extern "C" {
#ifdef HAVE_MA86
#include "hsl_ma86_wrapper.h"
#endif
#ifdef HAVE_QDLDL
#include "qdldl.h"
#endif
}

enum DecomposerStatus {
  kDecomposerStatusMin = 0,
  kDecomposerStatusOk = kDecomposerStatusMin,
  kDecomposerStatusErrorOom,
  kDecomposerStatusErrorFactorize,
  kDecomposerStatusErrorSolve,
  kDecomposerStatusErrorClear,
  kDecomposerStatusMax = kDecomposerStatusErrorClear
};

const double ok_theta_relative_tolerance = 1e-3;

// set default solver
// 1: ssids
// 2: ma86
// 3: qdldl
// 4: cholmod
// 5 : highs
const int default_solver = 2;

struct SsidsData {
  void* akeep{nullptr};
  void* fkeep{nullptr};
#ifdef HAVE_SPRAL
  struct spral_ssids_options options;
  struct spral_ssids_inform inform;
#endif
  int clear();
};

struct MA86Data {
  void* keep;
#ifdef HAVE_MA86
  ma86_control_d control;
  ma86_info_d info;
#endif
  std::vector<int> order;
  void clear();
};

struct QDLDLData {
#ifdef HAVE_QDLDL
  QDLDL_int Ln;
  QDLDL_int* Lp;
  QDLDL_int* Li;
  QDLDL_float* Lx;
  QDLDL_float* D;
  QDLDL_float* Dinv;

  QDLDL_int* etree;
  QDLDL_int* Lnz;
  QDLDL_int sumLnz;

  QDLDL_int* iwork;
  QDLDL_bool* bwork;
  QDLDL_float* fwork;

  QDLDL_float* x;  // Data for results of A\b
#endif
  void clear();
};

struct CholmodData {
#ifdef HAVE_CHOLMOD
  cholmod_common c;
  cholmod_triplet* T;
  cholmod_sparse* a;
  cholmod_factor* L;
  cholmod_dense* x;
  cholmod_dense* b;
#endif
  void clear();
};

struct HighsData {
  std::vector<HighsInt> basic_index;
  HFactor factor;
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
  HighsData highs_data;
  void clear(const int solver_type = default_solver);
};

void chooseDenseColumns(const HighsSparseMatrix& highs_a,
                        const std::vector<double>& theta,
                        const int option_max_dense_col,
                        const double option_dense_col_tolerance,
                        std::vector<int>& dense_col,
                        ExperimentData& experiment_data,
                        const bool quiet = true);

int augmentedInvert(const HighsSparseMatrix& highs_a,
                    const std::vector<double>& theta, IpmInvert& invert,
                    ExperimentData& experiment_data,
                    const int solver_type = default_solver);

void augmentedSolve(const HighsSparseMatrix& highs_a,
                    const std::vector<double>& theta,
                    const std::vector<double>& rhs_x,
                    const std::vector<double>& rhs_y,
                    std::vector<double>& lhs_x, std::vector<double>& lhs_y,
                    IpmInvert& invert, ExperimentData& experiment_data,
                    const int solver_type = default_solver);

double augmentedCondition(const HighsSparseMatrix& matrix,
                          const std::vector<double>& theta, IpmInvert& invert,
                          const int solver_type = default_solver);

int newtonInvert(const HighsSparseMatrix& highs_a,
                 const std::vector<double>& theta, IpmInvert& invert,
                 const int option_max_dense_col,
                 const double option_dense_col_tolerance,
                 ExperimentData& experiment_data, const bool quiet = true,
                 const int solver_type = default_solver);

int newtonSolve(const HighsSparseMatrix& highs_a,
                const std::vector<double>& theta,
                const std::vector<double>& rhs, std::vector<double>& lhs,
                IpmInvert& invert, ExperimentData& experiment_data,
                const int& solver_type = default_solver);

double newtonCondition(const HighsSparseMatrix& matrix,
                       const std::vector<double>& theta, IpmInvert& invert,
                       const int& solver_type = default_solver);

bool increasingIndex(const HighsSparseMatrix& matrix);
void productAThetaAT(const HighsSparseMatrix& matrix,
                     const std::vector<double>& theta,
                     const std::vector<double>& x, std::vector<double>& result);

int computeAThetaAT(const HighsSparseMatrix& matrix,
                    const std::vector<double>& theta, HighsSparseMatrix& AAT,
                    const int max_num_nz = 100000000
                    // Cant exceed kHighsIInf = 2,147,483,647,
                    // otherwise start_ values may overflow. Even
                    // 100,000,000 is probably too large, unless the
                    // matrix is near-full, since fill-in will
                    // overflow pointers
);

int analyseScaledRowNorms(const HighsSparseMatrix& matrix,
                          const std::vector<double>& theta,
                          const bool quiet = true);

int gepp(const std::vector<std::vector<double>>& matrix,
         const std::vector<double>& rhs, std::vector<double>& solution);

int callSsidsAugmentedFactor(const HighsSparseMatrix& matrix,
                             const std::vector<double>& theta,
                             SsidsData& ssids_data,
                             ExperimentData& experiment_data);
int callSsidsNewtonFactor(const HighsSparseMatrix& AThetaAT,
                          SsidsData& ssids_data,
                          ExperimentData& experiment_data);
void callSsidsSolve(const int system_size, const int num_rhs, double* rhs,
                    SsidsData& ssids_data);

int callMA86AugmentedFactor(const HighsSparseMatrix& matrix,
                            const std::vector<double>& theta,
                            MA86Data& ma86_data,
                            ExperimentData& experiment_data);

int callMA86NewtonFactor(const HighsSparseMatrix& AThetaAT, MA86Data& ma86_data,
                         ExperimentData& experiment_data);

void callMA86Solve(const int system_size, const int num_rhs, double* rhs,
                   MA86Data& ma86_data);

int callQDLDLAugmentedFactor(const HighsSparseMatrix& matrix,
                             const std::vector<double>& theta,
                             QDLDLData& qdldl_data,
                             ExperimentData& experiment_data);

int callQDLDLNewtonFactor(const HighsSparseMatrix& AThetaAT,
                          QDLDLData& qdldl_data,
                          ExperimentData& experiment_data);

void callQDLDLSolve(const int system_size, const int num_rhs, double* rhs,
                    QDLDLData& qdldl_data);

int callCholmodNewtonFactor(const HighsSparseMatrix& AThetaAT,
                            CholmodData& cholmod_data,
                            ExperimentData& experiment_data,
                            const int num_thods = 3);
void callCholmodSolve(const int system_size, const int num_rhs, double* rhs,
                      CholmodData& cholmod_data);

int callHighsAugmentedFactor(const HighsSparseMatrix& matrix,
                             const std::vector<double>& theta,
                             HighsData& highs_data,
                             ExperimentData& experiment_data);

int callHighsNewtonFactor(const HighsSparseMatrix& AThetaAT,
                          HighsData& highs_data,
                          ExperimentData& experiment_data);

void callHighsSolve(std::vector<std::vector<double>>& rhs,
                    HighsData& highs_data);

void callHighsSolve(std::vector<double>& rhs, HighsData& highs_data);

// ------------------------------------------
// Metis invert and solve
// ------------------------------------------
int blockInvert(const HighsSparseMatrix& block, IpmInvert& invert,
                ExperimentData& experiment_data,
                const int solver_type = default_solver);

int blockSolve(double* rhs, int num_rhs, IpmInvert& invert,
               ExperimentData& experiment_data,
               const int& solver_type = default_solver);

int diagonalForwardSolve(double* rhs, int nrhs, IpmInvert& invert,
                         ExperimentData& experiment_data,
                         double* DiagForwardSolvedRhs,
                         const int& solver_type = default_solver);

#endif
