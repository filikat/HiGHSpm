#include "IPM_aux.h"

#include <fstream>
#include <iostream>

#include "VectorOperations.h"

Residuals::Residuals(int m, int n)
    : res1(m, 0.0),
      res2(n, 0.0),
      res3(n, 0.0),
      res4(n, 0.0),
      res5(n, 0.0),
      res6(n, 0.0) {}

void Residuals::print(int iter) const {
  std::ofstream out_file;

  char str[50];

  snprintf(str, 50, "res1_%d.txt", iter);
  out_file.open(str);
  for (double d : res1) {
    out_file << d << '\n';
  }
  out_file.close();

  snprintf(str, 50, "res2_%d.txt", iter);
  out_file.open(str);
  for (double d : res2) {
    out_file << d << '\n';
  }
  out_file.close();

  snprintf(str, 50, "res3_%d.txt", iter);
  out_file.open(str);
  for (double d : res3) {
    out_file << d << '\n';
  }
  out_file.close();

  snprintf(str, 50, "res4_%d.txt", iter);
  out_file.open(str);
  for (double d : res4) {
    out_file << d << '\n';
  }
  out_file.close();

  snprintf(str, 50, "res5_%d.txt", iter);
  out_file.open(str);
  for (double d : res5) {
    out_file << d << '\n';
  }
  out_file.close();

  snprintf(str, 50, "res6_%d.txt", iter);
  out_file.open(str);
  for (double d : res6) {
    out_file << d << '\n';
  }
  out_file.close();
}

bool Residuals::isNaN() const {
  if (isnan(res1) || isnan(res2) || isnan(res3) || isnan(res4) || isnan(res5) ||
      isnan(res6))
    return true;
  return false;
}

bool Residuals::isInf() const {
  if (isinf(res1) || isinf(res2) || isinf(res3) || isinf(res4) || isinf(res5) ||
      isinf(res6))
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

bool Iterate::isInf() const {
  if (isinf(x) || isinf(xl) || isinf(xu) || isinf(y) || isinf(zl) || isinf(zu))
    return true;
  return false;
}

void Iterate::print(int iter) const {
  std::ofstream out_file;

  char str[50];

  snprintf(str, 50, "it_x_%d.txt", iter);
  out_file.open(str);
  for (double d : x) {
    out_file << d << '\n';
  }
  out_file.close();

  snprintf(str, 50, "it_y_%d.txt", iter);
  out_file.open(str);
  for (double d : y) {
    out_file << d << '\n';
  }
  out_file.close();

  snprintf(str, 50, "it_xl_%d.txt", iter);
  out_file.open(str);
  for (double d : xl) {
    out_file << d << '\n';
  }
  out_file.close();

  snprintf(str, 50, "it_xu_%d.txt", iter);
  out_file.open(str);
  for (double d : xu) {
    out_file << d << '\n';
  }
  out_file.close();

  snprintf(str, 50, "it_zl_%d.txt", iter);
  out_file.open(str);
  for (double d : zl) {
    out_file << d << '\n';
  }
  out_file.close();

  snprintf(str, 50, "it_zu_%d.txt", iter);
  out_file.open(str);
  for (double d : zu) {
    out_file << d << '\n';
  }
  out_file.close();
}

NewtonDir::NewtonDir(int m, int n)
    : x(n, 0.0), y(m, 0.0), xl(n, 0.0), xu(n, 0.0), zl(n, 0.0), zu(n, 0.0) {}

bool NewtonDir::isNaN() const {
  if (isnan(x) || isnan(xl) || isnan(xu) || isnan(y) || isnan(zl) || isnan(zu))
    return true;
  return false;
}

bool NewtonDir::isInf() const {
  if (isinf(x) || isinf(xl) || isinf(xu) || isinf(y) || isinf(zl) || isinf(zu))
    return true;
  return false;
}

void NewtonDir::print(int iter) const {
  std::ofstream out_file;

  char str[50];

  snprintf(str, 50, "D_x_%d.txt", iter);
  out_file.open(str);
  for (double d : x) {
    out_file << d << '\n';
  }
  out_file.close();

  snprintf(str, 50, "D_y_%d.txt", iter);
  out_file.open(str);
  for (double d : y) {
    out_file << d << '\n';
  }
  out_file.close();

  snprintf(str, 50, "D_xl_%d.txt", iter);
  out_file.open(str);
  for (double d : xl) {
    out_file << d << '\n';
  }
  out_file.close();

  snprintf(str, 50, "D_xu_%d.txt", iter);
  out_file.open(str);
  for (double d : xu) {
    out_file << d << '\n';
  }
  out_file.close();

  snprintf(str, 50, "D_zl_%d.txt", iter);
  out_file.open(str);
  for (double d : zl) {
    out_file << d << '\n';
  }
  out_file.close();

  snprintf(str, 50, "D_zu_%d.txt", iter);
  out_file.open(str);
  for (double d : zu) {
    out_file << d << '\n';
  }
  out_file.close();
}
