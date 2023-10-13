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

// [ * *         + + ]
// [ * *         + + ]
// [     * *     & & ]
// [     * *     & & ]
// [         * * + + ]
// [         * * + + ]
// [ + + & & + + % % ]
// [ + + & & + + % % ]
//
// *    : diagonal blocks
// +, & : linking blocks
// %    : schur block (it will become the Schur complement)

class Metis_caller {
  // adjacency matrix of the graph
  HighsSparseMatrix M;

  // constraint matrix
  const HighsSparseMatrix* A;

  // type of M (augmented or normal equations)
  MetisPartitionType type;

  // number of parts in the partition
  int nparts;

  // number of vertices and edges in the graph
  int nvertex;
  int nedges;

  // partition of the graph.
  // partition[i] is in [0, nparts - 1] when Metis returns.
  // partition[i] is in [0, nparts] after the permutation is computed, because
  //   the linking block is counted as another part.
  std::vector<int> partition;

  // permutation and inverse permutation for matrix M
  std::vector<int> permutation;
  std::vector<int> perminv;

  // Contains the size of the blocks, after getPermutation returns.
  // blockSize[i], 0 <= i < nparts: size of the diagonal blocks
  // blockSize[nparts]: size of the linking block
  std::vector<int> blockSize;

  // Contains the blocks, after getBlocks returns.
  // Blocks[2 * i] contains the i-th diagonal block
  // Blocks[2 * i + 1] contains the i-th linking block
  // Blocks[2 * nparts] contains the Schur block (i.e. the last diagonal
  //   block that will become the Schur compl)
  std::vector<HighsSparseMatrix> Blocks;

  // Number of nonzeros of each of the diagonal and linking blocks, after
  // getNonzeros returns.
  // nzCount[2 * i] contains the nonzero count of the i-th diag block
  // nzCount[2 * i + 1] contains the nonzero count of the i-th linking block
  // nzCount[2 * nparts] contains the nonzero count of the Schur block
  // nzCount[2 * nparts + 1] should not be used for anything
  std::vector<int> nzCount;

  // if debug is true, print out multiple files for debug, using out_file
  bool debug = false;
  std::ofstream out_file;

  // store the factorization and data of the diagonal blocks
  std::vector<IpmInvert> invertData;
  std::vector<ExperimentData> expData;

 public:
  // Constructor
  // Set up matrix M with either augmented system or normal equations, depending
  // on the value of type.
  Metis_caller(const HighsSparseMatrix& input_A, MetisPartitionType input_type,
               int input_nparts);

  // Call Metis and produce the partition of the graph
  void getPartition();

  // From the partition, obtain the permutation to use for matrix M.
  // Two approached are tried to obtain a vertex cover of the edge cut:
  // - find a maximal matching of the edge cut (10 times)
  // - "greedy" heuristic: for each edge in the cut, include in the vertex cover
  //    the node with highest original numbering (when possible)
  // The one which yields the smallest Schur complement is used.
  void getPermutation();

  // Extracts the diagonal and linking blocks of matrix M permuted.
  // For augm system, diag1 and diag2 are the diagonals of 1,1 and 2,2 blocks.
  // For normal equations, diag1 is theta, diag2 is ignored.
  // At the first call, Blocks are computed. At subsequent calls:
  // - Blocks are recomputed entirely for normal equations
  // - only diagonal elements are updated for augm system
  void getBlocks(const std::vector<double>& diag1,
                 const std::vector<double>& diag2);

  // factorize and solve with the diagonal blocks
  // (temporary stuff for testing, schur complement still missing)
  void factor();
  void solve();

  void setDebug(bool db = true) { debug = db; }
  const HighsSparseMatrix& accessBlock(int i) const { return Blocks[i]; }

 private:
  // Computes the number of nonzeros of each of the diagonal and linking blocks.
  void getNonzeros();

  // Update diagonal elements of the augmented system blocks
  void updateDiag(const std::vector<double>& diag1,
                  const std::vector<double>& diag2);

  // print for debug
  template <typename T>
  void debug_print(const std::vector<T>& vec, const std::string& filename);
  void debug_print(const HighsSparseMatrix& mat, const std::string& filename);
};

template <typename T>
void Metis_caller::debug_print(const std::vector<T>& vec,
                               const std::string& filename) {
  out_file.open(filename);
  for (auto a : vec) {
    out_file << a << '\n';
  }
  out_file.close();
}

void denseMatrixPlusVector(std::vector<std::vector<double>>& mat,
                           const std::vector<double>& vec, int col);

#endif