#include "Symbolic.h"

#include <iostream>

#include "FactorHiGHSSettings.h"
#include "auxiliary/HpmLog.h"

namespace highspm {

Symbolic::Symbolic() : block_size_{kBlockSize} {}

void Symbolic::setParallel(bool par_tree, bool par_node) {
  parallel_tree_ = par_tree;
  parallel_node_ = par_node;
}

int64_t Symbolic::nz() const { return nz_; }
double Symbolic::flops() const { return flops_; }
double Symbolic::spops() const { return spops_; }
double Symbolic::critops() const { return critops_; }
Int Symbolic::blockSize() const { return block_size_; }
Int Symbolic::size() const { return n_; }
Int Symbolic::sn() const { return sn_; }
Int Symbolic::rows(Int i) const { return rows_[i]; }
Int Symbolic::ptr(Int i) const { return ptr_[i]; }
Int Symbolic::snStart(Int i) const { return sn_start_[i]; }
Int Symbolic::snParent(Int i) const { return sn_parent_[i]; }
Int Symbolic::relindCols(Int i) const { return relind_cols_[i]; }
Int Symbolic::relindClique(Int i, Int j) const { return relind_clique_[i][j]; }
Int Symbolic::consecutiveSums(Int i, Int j) const {
  return consecutive_sums_[i][j];
}
Int Symbolic::cliqueBlockStart(Int sn, Int bl) const {
  return clique_block_start_[sn][bl];
}
Int Symbolic::cliqueSize(Int sn) const {
  return clique_block_start_[sn].back();
}
bool Symbolic::parTree() const { return parallel_tree_; }
bool Symbolic::parNode() const { return parallel_node_; }

const std::vector<Int>& Symbolic::ptr() const { return ptr_; }
const std::vector<Int>& Symbolic::iperm() const { return iperm_; }
const std::vector<Int>& Symbolic::snParent() const { return sn_parent_; }
const std::vector<Int>& Symbolic::snStart() const { return sn_start_; }
const std::vector<Int>& Symbolic::pivotSign() const { return pivot_sign_; }

void printMemory(double mem) {
  if (mem < 1024)
    Log::printf("serial memory   : %.1f B\n", mem);
  else if (mem < 1024 * 1024)
    Log::printf("serial memory   : %.1f KB\n", mem / 1024);
  else if (mem < 1024 * 1024 * 1024)
    Log::printf("serial memory   : %.1f MB\n", mem / 1024 / 1024);
  else
    Log::printf("serial memory   : %.1f GB\n", mem / 1024 / 1024 / 1024);
}
void Symbolic::print(bool verbose) const {
  Log::printf("\nSymbolic factorisation:\n");
  Log::printf("size            : %.2e\n", (double)n_);
  Log::printf("nnz             : %.2e\n", (double)nz_);
  Log::printf("fill-in         : %.2f\n", fillin_);
  printMemory(serial_storage_);
  Log::printf("flops           : %.1e\n", flops_);

  if (verbose) {
    Log::printf("sparse ops      : %.1e\n", spops_);
    Log::printf("critical ops    : %.1e\n", critops_);
    Log::printf("max tree speedup: %.2f\n", flops_ / critops_);
    Log::printf("artificial nz   : %.1e (%.1f%%)\n", (double)artificial_nz_,
                (double)artificial_nz_ / nz_ * 100);
    Log::printf("artificial ops  : %.1e (%.1f%%)\n", artificial_ops_,
                artificial_ops_ / flops_ * 100);
    Log::printf("largest front   : %5d\n", largest_front_);
    Log::printf("largest sn      : %5d\n", largest_sn_);
    Log::printf("supernodes      : %5d\n", sn_);
    Log::printf("sn size <=   1  : %5d\n", sn_size_1_);
    Log::printf("sn size <=  10  : %5d\n", sn_size_10_);
    Log::printf("sn size <= 100  : %5d\n", sn_size_100_);
    Log::printf("sn avg size     : %5.1f\n", (double)n_ / sn_);
  }

  Log::printf("\n");
}

}  // namespace highspm