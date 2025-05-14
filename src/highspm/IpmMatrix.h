#ifndef HIGHSPM_IPM_MATRIX_H
#define HIGHSPM_IPM_MATRIX_H

#include "auxiliary/IntConfig.h"
#include "auxiliary/KrylovMethods.h"
#include "util/HighsSparseMatrix.h"

namespace highspm {

class IpmMatrix : public AbstractMatrix {
  const HighsSparseMatrix* A_;
  const double* scaling_;
  const Int n_;
  bool use_as_;

 public:
  IpmMatrix(const HighsSparseMatrix& A, const std::vector<double>& scaling,
            bool use_as);

  void apply(std::vector<double>& x) const override;
  void apply(std::vector<HighsCDouble>& x) const override;

  void applyNE(std::vector<double>& x) const;
  void applyNE(std::vector<HighsCDouble>& x) const;
  void applyAS(std::vector<double>& x) const;
  void applyAS(std::vector<HighsCDouble>& x) const;
};

}  // namespace highspm

#endif