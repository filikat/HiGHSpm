#ifndef FACTORHIGHS_FACTORISE_H
#define FACTORHIGHS_FACTORISE_H

#include <cmath>

#include "Numeric.h"
#include "Symbolic.h"

class Factorise {
 public:
  // matrix to factorise
  std::vector<int> rowsA_{};
  std::vector<int> ptrA_{};
  std::vector<double> valA_{};
  int n_{};
  int nzA_{};

  // symbolic factorisation
  const Symbolic& S_;

  // children in supernodal elimination tree
  std::vector<int> first_child_{};
  std::vector<int> next_child_{};

  // reverse linked lists of chidlren
  std::vector<int> first_child_reverse_{};
  std::vector<int> next_child_reverse_{};

  // generated elements, aka Schur complements.
  std::vector<std::vector<double>> schur_contribution_{};

  // columns of L, stored as dense supernodes
  std::vector<std::vector<double>> sn_columns_{};

  // swaps of columns for each supernode, ordered locally within a block
  std::vector<std::vector<int>> swaps_{};

  // Information about 2x2 pivots.
  // If pivot_2x2[sn][i] == 0, 1x1 pivot was used.
  // If pivot_2x2[sn][i] != 0, 2x2 pivot was used and pivot_2x2[sn][i] stores
  //  the off-diagonal pivot entry (of the 2x2 inverse).
  std::vector<std::vector<double>> pivot_2x2_{};

  // largest diagonal element in the original matrix
  double max_diag_{};
  double min_diag_{};
  double A_norm1_{};

  // regularization
  std::vector<double> total_reg_{};

  // flag to stop computation
  bool flag_stop_ = false;

 public:
  void permute(const std::vector<int>& iperm);
  void processSupernode(int sn);

 public:
  Factorise(const Symbolic& S, const std::vector<int>& rowsA,
            const std::vector<int>& ptrA, const std::vector<double>& valA);

  bool run(Numeric& num);
};

#endif