#include "IPM_aux.h"
#include "VectorOperations.h"
#include <iostream>
#include <cassert>

Residuals::Residuals(int m, int n)
    : res1(m, 0.0), res2(n, 0.0), res3(n, 0.0), res4(n, 0.0), res5(n, 0.0),
      res6(n, 0.0) {}

void Residuals::Print() const {
  std::cout << "\nResidual 1\n";
  for (double d : res1) {
    std::cout << d << ' ';
  }
  std::cout << "\nResidual 2\n";
  for (double d : res2) {
    std::cout << d << ' ';
  }
  std::cout << "\nResidual 3\n";
  for (double d : res3) {
    std::cout << d << ' ';
  }
  std::cout << "\nResidual 4\n";
  for (double d : res4) {
    std::cout << d << ' ';
  }
  std::cout << "\nResidual 5\n";
  for (double d : res5) {
    std::cout << d << ' ';
  }
  std::cout << "\nResidual 6\n";
  for (double d : res6) {
    std::cout << d << ' ';
  }
  std::cout << '\n';
}

bool Residuals::isNaN() const {
  if (isnan(res1) || isnan(res2) || isnan(res3) || isnan(res4) || isnan(res5) ||
      isnan(res6))
    return true;
  return false;
}

Iterate::Iterate(int m, int n)
    : x(n, 0.0), y(m, 0.0), xl(n, 1.0), xu(n, 1.0), zl(n, 1.0), zu(n, 1.0) {}

bool Iterate::isNaN() const {
  if (isnan(x) || isnan(xl) || isnan(xu) || isnan(y) || isnan(zl) || isnan(zu))
    return true;
  return false;
}

NewtonDir::NewtonDir(int m, int n)
    : x(n, 0.0), y(m, 0.0), xl(n, 0.0), xu(n, 0.0), zl(n, 0.0), zu(n, 0.0) {}

bool NewtonDir::isNaN() const {
  if (isnan(x) || isnan(xl) || isnan(xu) || isnan(y) || isnan(zl) || isnan(zu))
    return true;
  return false;
}

std::string decomposerSource(int decomposer_source) {
  assert(decomposer_source >= kDecomposerSourceMin &&
	 decomposer_source <= kDecomposerSourceMax);
  if (decomposer_source == kDecomposerSourceSsids) {
    return "Ssids";
  } else if (decomposer_source == kDecomposerSourceMa86) {
    return "MA86";
  } else if (decomposer_source == kDecomposerSourceQdldl) {
    return "QDLDL";
  } else if (decomposer_source == kDecomposerSourceCholmod) {
    return "Cholmod";
  } else if (decomposer_source == kDecomposerSourceHighs) {
    return "HiGHS";
  } else  {
    return "Unknown";
  } 
}

