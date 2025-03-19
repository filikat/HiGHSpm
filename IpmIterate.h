#ifndef IPM_ITERATE_H
#define IPM_ITERATE_H

#include <vector>

#include "IpmModel.h"

// Holds the Newton direction Delta(x,y,xl,xu,zl,zu)
struct NewtonDir {
  std::vector<double> x{};
  std::vector<double> y{};
  std::vector<double> xl{};
  std::vector<double> xu{};
  std::vector<double> zl{};
  std::vector<double> zu{};

  NewtonDir(int m, int n);
};

struct IpmIterate {
  // lp model
  const IpmModel& model_;

  // ipm point
  std::vector<double> x_, xl_, xu_, y_, zl_, zu_;

  // residuals
  std::vector<double> res1_, res2_, res3_, res4_, res5_, res6_;

  // Newton direction
  NewtonDir delta_;
  std::vector<double>&dx_, &dxl_, &dxu_, &dy_, &dzl_, &dzu_;

  // indicators
  double pobj_, dobj_, pinf_, dinf_, pdgap_;

  double mu_;
  std::vector<double> scaling_;

  // ===================================================================================
  // Functions to construct, clear and check for nan or inf
  // ===================================================================================
  IpmIterate(const IpmModel& model);

  // clear existing data
  void clearIter();
  void clearRes();
  void clearDir();

  // check if any component is nan or infinite
  bool isNan() const;
  bool isInf() const;
  bool isResNan() const;
  bool isResInf() const;
  bool isDirNan() const;
  bool isDirInf() const;

  // ===================================================================================
  // Compute:
  //  mu = \sum xl(i) * zl(i) + \sum xu(j) * zu(j)
  // for variables: i with a finite lower bound
  //                j with a finite upper bound
  // ===================================================================================
  void mu();

  // ===================================================================================
  // Compute diagonal scaling Theta^{-1}
  //  Theta^{-1}_{ii} = zl(i) / xl(i) + zu(i) / xu(i)
  // Theta^{-1} only considers the terms above if the corresponding upper/lower
  // bound is finite.
  // ===================================================================================
  void scaling();

  // ===================================================================================
  // Compute convergence indicators
  // - primal and dual objectives, and primal-dual relative gap
  // - primal and dual infeasibilities (either scaled or unscaled)
  // - complementairy products
  // ===================================================================================
  void indicators();

  // objectives
  void primalObj();
  void dualObj();
  void pdGap();

  // infeasibilities
  void primalInfeas();
  void dualInfeas();
  void primalInfeasUnscaled();
  void dualInfeasUnscaled();
  double infeasAfterDropping() const;

  // complementarity products
  void products();

  // ===================================================================================
  // Compute:
  //  res1 = rhs - A * x
  //  res2 = lower - x + xl
  //  res3 = upper - x - xu
  //  res4 = c - A^T * y - zl + zu
  // Components of residuals 2,3 are set to zero if the corresponding
  // upper/lower bound is not finite.
  // ===================================================================================
  void residual1234();

  // ===================================================================================
  // Compute:
  //  res5 = sigma * mu * e - Xl * Zl * e
  //  res6 = sigma * mu * e - Xu * Zu * e
  // Components of residuals 5,6 are set to zero if the corresponding
  // upper/lower bound is not finite.
  // ===================================================================================
  void residual56(double sigma);

  // ===================================================================================
  // Compute:
  //  res7 = res4 - Xl^{-1} * (res5 + Zl * res2) + Xu^{-1} * (res6 - Zu * res3)
  // (the computation of res7 takes into account only the components for which
  // the correspoding upper/lower bounds are finite)
  // ===================================================================================
  std::vector<double> residual7() const;

  // ===================================================================================
  // Compute:
  //  res8 = res1 + A * Theta * res7
  // ===================================================================================
  std::vector<double> residual8(const std::vector<double>& res7) const;

  // ===================================================================================
  // Extract solution to be returned to user:
  // - remove extra slacks from x, xl, xu, zl, zu
  // - adjust sign of y for inequality constraints
  // - compute and adjust sign of slacks
  // ===================================================================================
  void extract(std::vector<double>& x, std::vector<double>& xl,
               std::vector<double>& xu, std::vector<double>& slack,
               std::vector<double>& y, std::vector<double>& zl,
               std::vector<double>& zu) const;

  // ===================================================================================
  // Extract complementary solution to be used for crossover with IPX:
  // - drop variables to obtain complementary (x,y,z)
  // - adjust y based on z-slacks
  // - compute slacks
  // - remove extra slacks from x, z
  // ===================================================================================
  void extract(std::vector<double>& x, std::vector<double>& slack,
               std::vector<double>& y, std::vector<double>& z) const;

  // ===================================================================================
  // Construct a complementary point (x,y,z), such that for each j, either xj is
  // at one of the bounds (lower or upper), or zj is zero.
  // ===================================================================================
  void dropToComplementarity(std::vector<double>& x, std::vector<double>& y,
                             std::vector<double>& z) const;
};

#endif