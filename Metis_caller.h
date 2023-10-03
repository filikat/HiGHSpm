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

// Call Metis
// A: constraint matrix
// type: MetisPartitionType for augmented system or normal equations
// nparts: number of parts in metis partition
// permutation: on return, contains the permutation
void getMetisPermutation(const HighsSparseMatrix& A, MetisPartitionType type,
                         int nparts, std::vector<int>& permutation,
                         std::vector<int>& blockSize);

// Extracts the diagonal and linking blocks of matrix M with permutation.
// On return, Blocks[2 * i] contains the i-th diagonal block
//            Blocks[2 * i + 1] contains the i-th linking block
//            Blocks[2 * nparts + 1] contains the Schur block
//            (i.e. the last diagonal block that will become the Schur compl)
void getBlocks(const std::vector<int>& permutation,
               const std::vector<int>& blockSize, MetisPartitionType type,
               const HighsSparseMatrix& M,
               std::vector<HighsSparseMatrix>& Blocks);

// Given the inverse permutation perminv and the block sizes, computes the
// number of nonzeros of each of the diagonal and linking blocks.
// On return, nzCount[2 * i] contains the nonzero count of the i-th diag block
//            nzCount[2 * i + 1] contains the nonzero count of the i-th linking
//            block
//            nzCount[2 * nparts + 1] contains the nonzero count of the Schur
//            block
//            nzCount[2 * nparts + 2] should not be used for anything
void getNonzeros(const HighsSparseMatrix& M, const std::vector<int>& perminv,
                 const std::vector<int>& blockSize, std::vector<int>& nzCount);

#endif