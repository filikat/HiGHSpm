#ifndef IPM_AUX_H
#define IPM_AUX_H

#include "lp_data/HConst.h"
#include "IPM_const.h"

#include <string>
#include <vector>

const std::string kModelFileCommandString = "model";
const std::string kPresolveCommandString = "presolve";
const std::string kDecomposerCommandString = "decomposer";
const std::string kNlaCommandString = "nla";
const std::string kPredcorCommandString = "predcor";
const std::string kMaxDenseColCommandString = "max_dense_col";
const std::string kDenseColToleranceCommandString = "density";

struct IpmOptions {
  std::string model_file = "";
  std::string model_name = "";
  std::string presolve = kHighsOnString;
  int decomposer_source = kOptionDecomposerSourceDefault;
  int nla = kOptionNlaDefault;
  std::string predcor = kHighsOnString;
  int max_dense_col = kOptionMaxDenseColDefault;
  double dense_col_tolerance = kOptionDenseColToleranceDefault;
  // Now the options set internally
  int iteration_limit = 100;
  double ipm_tolerance = 1e-6;
};

std::string decomposerSource(int decomposer_source);
std::string systemSolved(int nla);

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
  Residuals(int m, int n);

  void Print() const;
  bool isNaN() const;

  friend class IPM_caller;
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

  friend class IPM_caller;
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

  friend class IPM_caller;
};

// =======================================================================
// OUTPUT
// =======================================================================
struct Output {

  Iterate It{};
  int iterations{};
  double primal_infeas{};
  double dual_infeas{};
  double mu{};
};

#endif
