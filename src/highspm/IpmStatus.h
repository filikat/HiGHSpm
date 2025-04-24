#ifndef HIGHSPM_IPM_STATUS_H
#define HIGHSPM_IPM_STATUS_H

namespace highspm {

enum IpmStatus {
  kIpmStatusError,
  kIpmStatusMaxIter,
  kIpmStatusNoProgress,
  kIpmStatusOptimal,
  kIpmStatusPDFeas,
  kIpmStatusBasic
};

enum LinearSolverStatus {
  kLinearSolverStatusOk = 0,
  kLinearSolverStatusErrorOom,
  kLinearSolverStatusErrorAnalyse,
  kLinearSolverStatusErrorFactorise,
  kLinearSolverStatusErrorSolve
};

}  // namespace highspm

#endif