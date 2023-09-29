#ifndef METIS_CALLER_H
#define METIS_CALLER_H

#include "Direct.h"
#include "GKlib.h"
#include "VertexCover.h"
#include "metis.h"
#include "util/HighsSparseMatrix.h"

enum MetisPartitionType {
  kMetisAugmented,
  kMetisNormalEq,
};

// call Metis
// A: constraint matrix
// type: 0 for augmented system
//       1 for normal equations
// nparts: number of parts in metis partition
// permutation: on return, contains the permutation
void getMetisPermutation(const HighsSparseMatrix& A, MetisPartitionType type,
                         int nparts, std::vector<int>& permutation,
                         std::vector<int>& blockSize);

void getBlocks(const std::vector<int>& permutation,
               const std::vector<int>& blockSize, MetisPartitionType type,
               const HighsSparseMatrix& M,
               std::vector<HighsSparseMatrix>& blocks);

#endif