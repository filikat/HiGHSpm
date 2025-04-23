#include "Symbolic.h"

#include <iostream>

#include "FactorHiGHSSettings.h"

namespace highspm {

Symbolic::Symbolic(FormatType format_type)
    : format_type_{format_type}, block_size_{kBlockSize} {}

void Symbolic::setParallel(bool par_tree, bool par_node) {
  parallel_tree_ = par_tree;
  parallel_node_ = par_node;
}

FormatType Symbolic::formatType() const { return format_type_; }
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

}  // namespace highspm