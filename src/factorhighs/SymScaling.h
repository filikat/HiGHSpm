#ifndef SYM_SCALING_H
#define SYM_SCALING_H

#include <cmath>
#include <vector>

// Scalings for symmetric matrices, provided as lower triangular.

void CurtisReidScalingSym(const std::vector<int>& ptr,
                          const std::vector<int>& rows,
                          const std::vector<double>& val,
                          std::vector<double>& colscale);

void RuizScalingSym(const std::vector<int>& ptr, const std::vector<int>& rows,
                    const std::vector<double>& val,
                    std::vector<double>& colscale);

void JacekScalingSym(const std::vector<int>& ptr, const std::vector<int>& rows,
                     const std::vector<double>& val,
                     std::vector<double>& colscale);

#endif