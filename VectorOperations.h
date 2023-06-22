#include <vector>

// =======================================================================
// COMPONENT-WISE VECTOR OPERATIONS
// =======================================================================

// v1[i] + alpha * v2[i]
void VectorAdd(std::vector<double> &v1, const std::vector<double> &v2,
               double alpha);

// v1[i] + alpha
void VectorAdd(std::vector<double> &v1, const double alpha);

// alpha * v1[i] * v2[i] + beta
void VectorMultiply(std::vector<double> &v1, const std::vector<double> &v2,
                    double alpha, double beta);

// v1[i] + alpha * v2[i] * v3[i]
void VectorAddMult(std::vector<double> &v1, const std::vector<double> &v2,
                   const std::vector<double> &v3, double alpha);

// v1[i] / v2[i]
void VectorDivide(std::vector<double> &v1, const std::vector<double> &v2);

// v1[i] * alpha
void VectorScale(std::vector<double> &v1, double alpha);

// =======================================================================

// scalar product
double DotProd(const std::vector<double> &v1, const std::vector<double> &v2);

// Euclidean norm
double Norm2(const std::vector<double> &x);
