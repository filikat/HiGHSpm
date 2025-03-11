#ifndef IPM_AUX_H
#define IPM_AUX_H

#include <string>
#include <vector>

#include "Ipm_const.h"
#include "util/HighsSparseMatrix.h"

enum DecomposerStatus {
  kDecomposerStatusMin = 0,
  kDecomposerStatusOk = kDecomposerStatusMin,
  kDecomposerStatusErrorOom,
  kDecomposerStatusErrorFactorise,
  kDecomposerStatusErrorSolve,
  kDecomposerStatusErrorClear,
  kDecomposerStatusMax = kDecomposerStatusErrorClear
};

enum IpmStatus {
  kIpmStatusOptimal,
  kIpmStatusError,
  kIpmStatusMaxIter,
  kIpmStatusNoProgress
};

int computeLowerAThetaAT(const HighsSparseMatrix& matrix,
                         const std::vector<double>& scaling,
                         HighsSparseMatrix& AAT,
                         const int max_num_nz = 100000000
                         // Cant exceed kHighsIInf = 2,147,483,647,
                         // otherwise start_ values may overflow. Even
                         // 100,000,000 is probably too large, unless the
                         // matrix is near-full, since fill-in will
                         // overflow pointers
);

void debug_print(std::string filestr, const std::vector<int>& data);
void debug_print(std::string filestr, const std::vector<double>& data);

#endif
