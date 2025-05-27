#ifndef HIGHSPM_IPM_STATUS_H
#define HIGHSPM_IPM_STATUS_H

namespace highspm {

enum IpmStatus {
  // With these statuses, the solver can proceed with refining, if requested
  kIpmStatusNotRun,
  kIpmStatusMaxIter,
  kIpmStatusNoProgress,

  // With these statuses, the solver should stop and not attempt refining
  kIpmStatusStop,
  kIpmStatusError,
  kIpmStatusTimeLimit,
  kIpmStatusUserInterrupt,
  kIpmStatusPrimalInfeasible,
  kIpmStatusDualInfeasible,

  // Solver is optimal
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