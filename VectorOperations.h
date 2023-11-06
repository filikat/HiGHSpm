#ifndef VERTEX_OPERATIONS_H
#define VERTEX_OPERATIONS_H

#include <vector>

// =======================================================================
// COMPONENT-WISE VECTOR OPERATIONS
// =======================================================================

// v1[i] + alpha * v2[i]
void VectorAdd(std::vector<double>& v1, const std::vector<double>& v2,
               double alpha = 1.0);

// v1[i] + alpha
void VectorAdd(std::vector<double>& v1, const double alpha);

// alpha * v1[i] * v2[i] + beta
void VectorMultiply(std::vector<double>& v1, const std::vector<double>& v2,
                    double alpha = 1.0, double beta = 0.0);

// v1[i] + alpha * v2[i] * v3[i]
void VectorAddMult(std::vector<double>& v1, const std::vector<double>& v2,
                   const std::vector<double>& v3, double alpha = 1.0);

// v1[i] / v2[i]
void VectorDivide(std::vector<double>& v1, const std::vector<double>& v2);

// v1[i] * alpha
void VectorScale(std::vector<double>& v1, double alpha);

// =======================================================================

// scalar product
double DotProd(const std::vector<double>& v1, const std::vector<double>& v2);

// Euclidean norm
double Norm2(const std::vector<double>& x);
double Norm2(const std::vector<double>& x0, const std::vector<double>& x1);

double infNorm(const std::vector<double>& x);

// Infinity norm of the difference of two vectors
double infNormDiff(const std::vector<double>& x, const std::vector<double>& y);

// check for NaN
bool isnan(const std::vector<double>& x);

// check for Inf
bool isinf(const std::vector<double>& x);

#endif