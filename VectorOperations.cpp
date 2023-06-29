#include "VectorOperations.h"
#include <cassert>
#include <cmath>

void VectorAdd(std::vector<double> &v1, const std::vector<double> &v2,
               double alpha) {
  for (int i = 0; i < v1.size(); ++i) {
    v1[i] += alpha * v2[i];
  }
}

void VectorAdd(std::vector<double> &v1, const double alpha) {
  for (int i = 0; i < v1.size(); ++i) {
    v1[i] += alpha;
  }
}

void VectorMultiply(std::vector<double> &v1, const std::vector<double> &v2,
                    double alpha, double beta) {
  for (int i = 0; i < v1.size(); ++i) {
    v1[i] = alpha * v1[i] * v2[i] + beta;
  }
}

void VectorAddMult(std::vector<double> &v1, const std::vector<double> &v2,
                   const std::vector<double> &v3, double alpha) {
  for (int i = 0; i < v1.size(); ++i) {
    v1[i] += alpha * v2[i] * v3[i];
  }
}

void VectorDivide(std::vector<double> &v1, const std::vector<double> &v2) {
  for (int i = 0; i < v1.size(); ++i) {
    v1[i] /= v2[i];
  }
}

void VectorScale(std::vector<double> &v1, double alpha) {
  for (int i = 0; i < v1.size(); ++i) {
    v1[i] *= alpha;
  }
}

double DotProd(const std::vector<double> &v1, const std::vector<double> &v2) {
  double result{};
  for (int i = 0; i < v1.size(); ++i) {
    result += v1[i] * v2[i];
  }
  return result;
}

double Norm2(const std::vector<double> &x) {
  double norm{};
  for (int i = 0; i < x.size(); ++i) {
    norm += (x[i] * x[i]);
  }
  return std::sqrt(norm);
}

double infNormDiff(const std::vector<double> &x, const std::vector<double> &y) {
  assert(x.size() == y.size());
  double inf_norm_diff = 0;
  for (int i = 0; i < int(x.size()); i++) {
    double diff = std::abs(x[i] - y[i]);
    inf_norm_diff = std::max(diff, inf_norm_diff);
  }
  return inf_norm_diff;
}

bool isnan(const std::vector<double> &x) {
  for (int i = 0; i < x.size(); ++i) {
    if (std::isnan(x[i]))
      return true;
  }
  return false;
}
