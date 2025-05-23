#ifndef HIGHSPM_HPM_INFO_H
#define HIGHSPM_HPM_INFO_H

#include "HpmOption.h"
#include "HpmStatus.h"
#include "auxiliary/IntConfig.h"
#include "ipm/ipx/info.h"

namespace highspm {

struct HpmInfo {
  // Size of problem, as seen by the solver
  Int m_solver, n_solver;

  // Size of original problem
  Int m_original, n_original;

  // Status of solver, see IpmStatus.h
  IpmStatus ipm_status = kIpmStatusNotRun;

  // Number of ipm iterations performed
  Int ipm_iter = 0;

  // True if ipx was invoked, whether to refine solution or for crossover
  bool ipx_used = false;

  // Info from ipx
  ipx::Info ipx_info;

  // Number of correctors used
  Int correctors;

  // Nla option used
  OptionNla option_nla;

  // Parallel option used
  OptionParallel option_par;

  // Total times to form matrix, factorise and solve linear systems
  double analyse_NE_time{};
  double analyse_AS_time{};
  double matrix_time{};
  double factor_time{};
  double solve_time{};

  // Counters
  Int factor_number{};
  Int solve_number{};

  // Information on dense columns
  Int num_dense_cols{};
  double max_col_density{};
};

}  // namespace highspm

#endif