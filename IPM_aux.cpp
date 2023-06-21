#include "IPM_aux.h"
#include <iostream>

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

Iterate::Iterate(int m, int n)
    : x(n, 0.0), y(m, 0.0), xl(n, 1.0), xu(n, 1.0), zl(n, 1.0), zu(n, 1.0) {}

NewtonDir::NewtonDir(int m, int n)
    : x(n, 0.0), y(m, 0.0), xl(n, 0.0), xu(n, 0.0), zl(n, 0.0), zu(n, 0.0) {}
