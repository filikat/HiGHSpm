#ifndef HIGHSPM_IPM_STATUS_H
#define HIGHSPM_IPM_STATUS_H

#include <map>
#include <string>

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

inline std::string statusString(IpmStatus status) {
  static const std::map<IpmStatus, std::string> status_map{
      {kIpmStatusNotRun, "not run"},
      {kIpmStatusMaxIter, "max iterations"},
      {kIpmStatusNoProgress, "no progress"},
      {kIpmStatusError, "internal error"},
      {kIpmStatusTimeLimit, "time limit"},
      {kIpmStatusUserInterrupt, "user interrupt"},
      {kIpmStatusPrimalInfeasible, "primal infeasible"},
      {kIpmStatusDualInfeasible, "dual infeasible"},
      {kIpmStatusPDFeas, "primal-dual feasible"},
      {kIpmStatusBasic, "crossover optimal"}};

  auto found = status_map.find(status);
  if (found != status_map.end()) return found->second;
  return "unknown";
}

}  // namespace highspm

#endif