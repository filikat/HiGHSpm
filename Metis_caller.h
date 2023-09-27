#ifndef METIS_CALLER_H
#define METIS_CALLER_H

#include "GKlib.h"
#include "metis.h"
#include "VertexCover.h"
#include "util/HighsSparseMatrix.h"
#include "Direct.h"


enum MetisPartitionType {
    kMetisAugmented,
    kMetisNormalEq,
};

// call Metis
// A: constraint matrix
// type: == 0 for augmented system
//       != 0 for normal equations
// nparts: number of parts in metis partition
// permutation: on return, contains the permutation
void getMetisPermutation(const HighsSparseMatrix& A, int type, int nparts,
                         std::vector<int>& permutation);

#endif