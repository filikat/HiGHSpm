#ifndef KRYLOV_METHODS_H
#define KRYLOV_METHODS_H

#include <vector>

#include "Numeric.h"
#include "util/HighsSparseMatrix.h"

// Abstract class for matrices inside of a Krylov method
class AbstractMatrix {
 public:
  virtual void apply(std::vector<double>& x) const = 0;
};

// Class to perform matrix-vector products with ipm matrix
class IpmMatrix : public AbstractMatrix {
  const HighsSparseMatrix* A_ = nullptr;
  const std::vector<double>* scaling_ = nullptr;
  bool use_as_ = false;

 public:
  void reset(const HighsSparseMatrix& A, const std::vector<double>& scaling,
             bool use_as = false);
  void apply(std::vector<double>& x) const override;
};

// Class to perform solves with the ipm factorization
class IpmFactor : public AbstractMatrix {
  const Numeric& N_;

 public:
  IpmFactor(const Numeric& N);
  void apply(std::vector<double>& x) const override;
};

// Class to apply diagonal preconditioner
class NeDiagPrec : public AbstractMatrix {
  std::vector<double> diag;

 public:
  void reset(const HighsSparseMatrix& A, const std::vector<double>& scaling);
  void apply(std::vector<double>& x) const override;
};

// =================
// Krylov solvers
// =================

int Gmres(const AbstractMatrix* M, const AbstractMatrix* P,
          const std::vector<double>& b, std::vector<double>& x, double tol,
          int maxit);

int Cg(const AbstractMatrix* M, const AbstractMatrix* P,
       const std::vector<double>& b, std::vector<double>& x, double tol,
       int maxit);

#endif