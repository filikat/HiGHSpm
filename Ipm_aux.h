#ifndef IPM_AUX_H
#define IPM_AUX_H

#include <string>
#include <vector>

#include "Ipm_const.h"
#include "util/HighsSparseMatrix.h"

enum DecomposerStatus {
  kDecomposerStatusMin = 0,
  kDecomposerStatusOk = kDecomposerStatusMin,
  kDecomposerStatusErrorOom,
  kDecomposerStatusErrorFactorise,
  kDecomposerStatusErrorSolve,
  kDecomposerStatusErrorClear,
  kDecomposerStatusMax = kDecomposerStatusErrorClear
};

// =======================================================================
// RESIDUALS
// =======================================================================
// Holds the six residuals of the IPM
// res1 = rhs - A * x
// res2 = lower - x + xl
// res3 = upper - x - xu
// res4 = obj - A^T * y - zl + zu
// res5 = sigma * mu * e - Xl * Zl * e
// res6 = sigma * mu * e - Xu * Zu * e
//
class Residuals {
  std::vector<double> res1{};
  std::vector<double> res2{};
  std::vector<double> res3{};
  std::vector<double> res4{};
  std::vector<double> res5{};
  std::vector<double> res6{};

 public:
  Residuals() = default;
  Residuals(int m, int n);

  void print(int iter = 0) const;
  bool isNaN() const;
  bool isInf() const;

  friend class Ipm;
};

// =======================================================================
// ITERATE
// =======================================================================
// Holds the iterate (x,y,xl,xu,zl,zu)
class Iterate {
 public:
  std::vector<double> x{};
  std::vector<double> y{};
  std::vector<double> xl{};
  std::vector<double> xu{};
  std::vector<double> zl{};
  std::vector<double> zu{};

  Iterate() = default;
  Iterate(int m, int n);

  bool isNaN() const;
  bool isInf() const;

  void print(int iter = 0) const;

  friend class Ipm;
};

// =======================================================================
// DIRECTION
// =======================================================================
// Holds the Newton direction Delta(x,y,xl,xu,zl,zu)
class NewtonDir {
  std::vector<double> x{};
  std::vector<double> y{};
  std::vector<double> xl{};
  std::vector<double> xu{};
  std::vector<double> zl{};
  std::vector<double> zu{};

 public:
  NewtonDir() = default;
  NewtonDir(int m, int n);

  bool isNaN() const;
  bool isInf() const;

  void print(int iter = 0) const;

  friend class Ipm;
};

// =======================================================================
// OUTPUT
// =======================================================================
struct Output {
  std::vector<double> x{};
  std::vector<double> xl{};
  std::vector<double> xu{};
  std::vector<double> slack{};
  std::vector<double> y{};
  std::vector<double> zl{};
  std::vector<double> zu{};
  int iterations{};
  double primal_infeas{};
  double dual_infeas{};
  double mu{};
  std::string status = "Error";
};

int computeAThetaAT(const HighsSparseMatrix& matrix,
                    const std::vector<double>& scaling, HighsSparseMatrix& AAT,
                    const int max_num_nz = 100000000
                    // Cant exceed kHighsIInf = 2,147,483,647,
                    // otherwise start_ values may overflow. Even
                    // 100,000,000 is probably too large, unless the
                    // matrix is near-full, since fill-in will
                    // overflow pointers
);

int computeLowerAThetaAT(const HighsSparseMatrix& matrix,
                         const std::vector<double>& scaling,
                         HighsSparseMatrix& AAT,
                         const int max_num_nz = 100000000
                         // Cant exceed kHighsIInf = 2,147,483,647,
                         // otherwise start_ values may overflow. Even
                         // 100,000,000 is probably too large, unless the
                         // matrix is near-full, since fill-in will
                         // overflow pointers
);

void debug_print(std::string filestr, const std::vector<int>& data);
void debug_print(std::string filestr, const std::vector<double>& data);

#endif
