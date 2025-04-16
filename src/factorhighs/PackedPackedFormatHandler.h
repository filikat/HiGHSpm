#ifndef FACTORHIGHS_PACKED_PACKED_FORMAT_H
#define FACTORHIGHS_PACKED_PACKED_FORMAT_H

#include "FormatHandler.h"

class PackedPackedFormatHandler : public FormatHandler {
  std::vector<Int> diag_start_;

  void initFrontal() override;
  void initClique() override;
  void assembleFrontal(Int i, Int j, double val) override;
  void assembleFrontalMultiple(Int num, const std::vector<double>& child,
                               Int nc, Int child_sn, Int row, Int col, Int i,
                               Int j) override;
  Int denseFactorise(double reg_thresh) override;
  void assembleClique(const std::vector<double>& child, Int nc,
                      Int child_sn) override;
  void extremeEntries() override;

 public:
  PackedPackedFormatHandler(const Symbolic& S, Int sn);
};

#endif