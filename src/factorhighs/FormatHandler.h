#ifndef FORMAT_HANDLER_H
#define FORMAT_HANDLER_H

#include <vector>

#include "Symbolic.h"

// Interface class to handle different formats of dense matrices during the
// factorise phase.
// Any implementation of a specific format needs to define:
// - initFrontal: to initialize the frontal matrix with the correct number of
//                elements; the entries should be set to zero.
// - initClique: to initialize the clique matrix with the correct number of
//                elements; the entries should be left uninitialized.
// - assembleFrontal: to set a specific entry of frontal (used to assemble the
//                original matrix)
// - assembleFrontalMultiple: to sum a given number of consecutive entries into
//                frontal (used to assemble the child supernodes)
// - assembleClique: to sum the contributions of a given child supernode into
//                clique.
// - denseFactorise: to perform the dense partial factorization of frontal,
//                storing the Schur complement in clique.

class FormatHandler {
 protected:
  // symbolic object
  const Symbolic* S_;

  // supernode being processed
  const int sn_{};

  // block size
  const int nb_{};

  // size of the supernode
  const int sn_size_{};

  // size of the front
  const int ldf_{};

  // size of the clique
  const int ldc_{};

  // local copies to be moved at the end
  std::vector<double> frontal_{};
  std::vector<double> clique_{};
  std::vector<double> local_reg_{};
  std::vector<int> swaps_{};
  std::vector<double> pivot_2x2_{};

 public:
  FormatHandler(const Symbolic& S, int sn);
  void terminate(std::vector<double>& frontal, std::vector<double>& clique,
                 std::vector<double>& total_reg, std::vector<int>& swaps,
                 std::vector<double>& pivot_2x2);

  // avoid copies
  FormatHandler(const FormatHandler&) = delete;
  FormatHandler& operator=(const FormatHandler&) = delete;

  // virtual destructor
  virtual ~FormatHandler() = default;

  // =================================================================
  // Pure virtual functions.
  // These need to be defined by any derived class.
  // =================================================================
  virtual void initFrontal() = 0;
  virtual void initClique() = 0;
  virtual void assembleFrontal(int i, int j, double val) = 0;
  virtual void assembleFrontalMultiple(int num,
                                       const std::vector<double>& child, int nc,
                                       int child_sn, int row, int col, int i,
                                       int j) = 0;
  virtual void assembleClique(const std::vector<double>& child, int nc,
                              int child_sn) = 0;
  virtual int denseFactorise(double reg_thresh) = 0;

  // =================================================================
  // Virtual functions.
  // These may be overridden by derived classes, if needed.
  // =================================================================
  virtual void extremeEntries() {}
};

const int extra_space = 10;

#endif