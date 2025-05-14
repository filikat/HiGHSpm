#ifndef HIGHSPM_KRYLOV_METHODS_H
#define HIGHSPM_KRYLOV_METHODS_H

#include <vector>

#include "auxiliary/IntConfig.h"
#include "util/HighsCDouble.h"

namespace highspm {

// Abstract class for matrices inside of a Krylov method
class AbstractMatrix {
 public:
  // virtual destructor
  virtual ~AbstractMatrix() = default;

  virtual void apply(std::vector<double>& x) const {};
  virtual void apply(std::vector<HighsCDouble>& x) const {};
};

// =================
// Krylov solvers
// =================

Int Gmres(const AbstractMatrix* M, const AbstractMatrix* P,
          const std::vector<double>& b, std::vector<double>& x, double tol,
          Int maxit);

Int Cg(const AbstractMatrix* M, const AbstractMatrix* P,
       const std::vector<double>& b, std::vector<double>& x, double tol,
       Int maxit);

}  // namespace highspm

#endif