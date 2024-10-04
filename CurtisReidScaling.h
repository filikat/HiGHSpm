#ifndef CURTIS_REID_SCALING_H
#define CURTIS_REID_SCALING_H

#include <cmath>
#include <vector>

#include "VectorOperations.h"

void CurtisReidScaling(const std::vector<int>& ptr,
                       const std::vector<int>& rows,
                       const std::vector<double>& val,
                       const std::vector<double> b, const std::vector<double> c,
                       double& objscale, double& rhsscale,
                       std::vector<double>& rowscale,
                       std::vector<double>& colscale);

#endif