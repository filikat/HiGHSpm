#ifndef HIGHSPM_IPM_INFO_H
#define HIGHSPM_IPM_INFO_H

#include "IpmStatus.h"
#include "auxiliary/IntConfig.h"
#include "ipm/ipx/info.h"

namespace highspm {

struct IpmInfo {
  // Size of problem, as seen by the solver
  Int m_solver, n_solver;

  // Size of original problem
  Int m_original, n_original;

  // Status of solver, see IpmStatus.h
  IpmStatus ipm_status = kIpmStatusNotRun;

  // Number of ipm iterations performed
  Int ipm_iter = 0;

  // True if ipx was invoked, wither to refine solution or for crossover
  bool ipx_used = false;

  // Info from ipx
  ipx::Info ipx_info;

  // Number of correctors used
  Int correctors;

  // Nla option used
  OptionNla option_nla;

  // Parallel option used
  OptionParallel option_par;
};

}  // namespace highspm

#endif