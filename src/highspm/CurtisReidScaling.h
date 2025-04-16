#ifndef HIGHSPM_CURTIS_REID_SCALING_H
#define HIGHSPM_CURTIS_REID_SCALING_H

#include <cmath>
#include <vector>

#include "auxiliary/IntConfig.h"
#include "auxiliary/VectorOperations.h"

namespace highspm {

void CurtisReidScaling(const std::vector<Int>& ptr,
                       const std::vector<Int>& rows,
                       const std::vector<double>& val, std::vector<Int>& rowexp,
                       std::vector<Int>& colexp);

}

#endif