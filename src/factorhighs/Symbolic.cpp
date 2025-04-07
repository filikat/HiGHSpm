#include "Symbolic.h"

#include <iostream>

#include "FactorHiGHSSettings.h"

Symbolic::Symbolic(FormatType format_type)
    : format_type_{format_type}, block_size_{kBlockSize} {}

FormatType Symbolic::formatType() const { return format_type_; }
double Symbolic::nz() const { return nz_; }
double Symbolic::flops() const { return flops_; }
double Symbolic::spops() const { return spops_; }
int Symbolic::blockSize() const { return block_size_; }
int Symbolic::size() const { return n_; }
int Symbolic::sn() const { return sn_; }
int Symbolic::rows(int i) const { return rows_[i]; }
int Symbolic::ptr(int i) const { return ptr_[i]; }
int Symbolic::snStart(int i) const { return sn_start_[i]; }
int Symbolic::snParent(int i) const { return sn_parent_[i]; }
int Symbolic::relindCols(int i) const { return relind_cols_[i]; }
int Symbolic::relindClique(int i, int j) const { return relind_clique_[i][j]; }
int Symbolic::consecutiveSums(int i, int j) const {
  return consecutive_sums_[i][j];
}
int Symbolic::cliqueBlockStart(int sn, int bl) const {
  return clique_block_start_[sn][bl];
}
int Symbolic::cliqueSize(int sn) const {
  return clique_block_start_[sn].back();
}

const std::vector<int>& Symbolic::ptr() const { return ptr_; }
const std::vector<int>& Symbolic::iperm() const { return iperm_; }
const std::vector<int>& Symbolic::snParent() const { return sn_parent_; }
const std::vector<int>& Symbolic::snStart() const { return sn_start_; }
const std::vector<int>& Symbolic::pivotSign() const { return pivot_sign_; }
