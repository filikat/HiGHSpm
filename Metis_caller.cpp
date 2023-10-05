#include "Metis_caller.h"

// extern "C" {
void setOptions(idx_t* options) { METIS_SetDefaultOptions(options); }
void callMetis(idx_t nvertex, idx_t nconstraints, idx_t* adj_ptr,
               idx_t* adj_lst, idx_t nparts, idx_t* options, idx_t* objval,
               idx_t* part) {
  idx_t status =
      METIS_PartGraphKway(&nvertex, &nconstraints, adj_ptr, adj_lst, NULL, NULL,
                          NULL, &nparts, NULL, NULL, options, objval, part);

  assert(status == METIS_OK);
}
//}

Metis_caller::Metis_caller(const HighsSparseMatrix& A,
                           MetisPartitionType input_type, int input_nparts) {
  nparts = input_nparts;
  type = input_type;

  // -----------------------------------------------------------
  // set up the augmented system
  // -----------------------------------------------------------
  if (type == kMetisAugmented) {
    nvertex = A.num_row_ + A.num_col_;
    nedges = A.numNz() * 2;

    // allocate space for augmented matrix
    M.num_col_ = nvertex;
    M.num_row_ = nvertex;
    M.start_.resize(nvertex + 1);
    M.index_.resize(nedges);
    M.value_.resize(nedges);

    // temporary A transpose
    HighsSparseMatrix A_t = A;
    A_t.ensureRowwise();

    // create pointers of augmented matrix
    for (int i = 0; i < A.num_col_ + 1; ++i) {
      M.start_[i] = A.start_[i];
    }
    int shift = A.num_col_;
    int ptr_shift = A.numNz();
    for (int i = 1; i < A.num_row_ + 1; ++i) {
      M.start_[i + shift] = A_t.start_[i] + ptr_shift;
    }

    // create adjacency list of augmented matrix
    for (int i = 0; i < A.numNz(); ++i) {
      M.index_[i] = A.index_[i] + A.num_col_;
      M.value_[i] = A.value_[i];
    }
    int adj_shift = A.numNz();
    for (int i = 0; i < A.numNz(); ++i) {
      M.index_[i + adj_shift] = A_t.index_[i];
      M.value_[i + adj_shift] = A_t.value_[i];
    }
  }
  // -----------------------------------------------------------
  // set up the normal equations
  // -----------------------------------------------------------
  else if (type == kMetisNormalEq) {
    std::vector<double> theta(A.num_col_, 1.0);
    computeAThetaAT(A, theta, M);
    nvertex = A.num_row_;
    nedges = M.numNz();
  } else {
    std::cerr << "Wront type of matrix for Metis parition\n";
    return;
  }

  // -----------------------------------------------------------
  // allocate space for later
  // -----------------------------------------------------------
  partition.resize(nvertex);
  perminv.resize(nvertex);
  Blocks.resize(2 * nparts + 1);
  nzCount.resize(2 * nparts + 2);

  // -----------------------------------------------------------
  // print for debug
  // -----------------------------------------------------------
  if (debug) {
    debug_print(A.start_, "debug_data/A_ptr.txt");
    debug_print(A.index_, "debug_data/A_adj.txt");
    debug_print(A.value_, "debug_data/A_val.txt");
  }
}

void Metis_caller::getMetisPartition() {
  // create space for outputs
  idx_t objval{};

  // initialize metis options
  idx_t options[METIS_NOPTIONS];
  setOptions(options);

  // call Metis to get the partition
  callMetis(nvertex, 1, M.start_.data(), M.index_.data(), nparts, options,
            &objval, partition.data());
}

void Metis_caller::getMetisPermutation() {
  // compute permutation with maximal matching
  std::vector<int> permutationMM(nvertex);
  std::vector<int> blockSizeMM(nparts + 1);
  vertexCoverMM(nvertex, nedges, nparts, partition, M.start_, M.index_,
                permutationMM, blockSizeMM);

  // compute permutation with greedy heuristic
  std::vector<int> permutationG(nvertex);
  std::vector<int> blockSizeG(nparts + 1);
  vertexCoverG(nvertex, nedges, nparts, partition, M.start_, M.index_,
               permutationG, blockSizeG);

  // print for debug
  if (debug) {
    debug_print(partition, "debug_data/partition.txt");
    debug_print(permutationMM, "debug_data/permMM.txt");
    debug_print(permutationG, "debug_data/permG.txt");
    debug_print(blockSizeMM, "debug_data/blockSizeMM.txt");
    debug_print(blockSizeG, "debug_data/blockSizeG.txt");
  }

  // select permutation with smallest Schur complement
  if (blockSizeMM.back() > blockSizeG.back()) {
    permutation = std::move(permutationG);
    blockSize = std::move(blockSizeG);
  } else {
    permutation = std::move(permutationMM);
    blockSize = std::move(blockSizeMM);
  }

  // update partition so that if node i is linking, partition[i] = nparts
  for (int i = 0; i < blockSize.back(); ++i) {
    partition[permutation[nvertex - 1 - i]] = nparts;
  }
}

void Metis_caller::getBlocks() {
  // get inverse permutation
  for (int i = 0; i < perminv.size(); ++i) {
    perminv[permutation[i]] = i;
  }

  int threshold = M.num_col_ - blockSize.back();

  // get number of nonzeros in blocks for preallocation
  getNonzeros();

  // file to write blocks
  std::ofstream out_file;

  // index of column to consider
  int colStart = 0;

  // go through the blocks
  for (int blockId = 0; blockId < nparts; ++blockId) {
    // indices to access the correct block in Blocks
    int diagBlockIndex = 2 * blockId;
    int linkBlockIndex = 2 * blockId + 1;

    // count nonzeros in diagonal block and linking block
    int current_nz_block{};
    int current_nz_link{};

    // allocate space for blocks
    Blocks[diagBlockIndex].start_.reserve(blockSize[blockId] + 1);
    Blocks[linkBlockIndex].start_.reserve(blockSize[blockId] + 1);
    Blocks[diagBlockIndex].index_.reserve(nzCount[2 * blockId]);
    Blocks[linkBlockIndex].index_.reserve(nzCount[2 * blockId + 1]);
    Blocks[diagBlockIndex].value_.reserve(nzCount[2 * blockId]);
    Blocks[linkBlockIndex].value_.reserve(nzCount[2 * blockId + 1]);

    // go through the columns in the order of permutation
    for (int i = colStart; i < colStart + blockSize[blockId]; ++i) {
      int col = permutation[i];

      if (type == kMetisAugmented) {
        // diagonal is not included in augmented system
        Blocks[diagBlockIndex].index_.push_back(i - colStart);
        Blocks[diagBlockIndex].value_.push_back(1.0);
        ++current_nz_block;
      }

      // go through the column
      for (int colEl = M.start_[col]; colEl < M.start_[col + 1]; ++colEl) {
        int permuted = perminv[M.index_[colEl]];

        // determine if current element belongs to diagonal block or linking
        // block
        if (permuted < threshold) {
          assert(permuted >= colStart &&
                 permuted < colStart + blockSize[blockId]);
          Blocks[diagBlockIndex].index_.push_back(permuted - colStart);
          Blocks[diagBlockIndex].value_.push_back(M.value_[colEl]);
          ++current_nz_block;
        } else {
          Blocks[linkBlockIndex].index_.push_back(permuted - threshold);
          Blocks[linkBlockIndex].value_.push_back(M.value_[colEl]);
          ++current_nz_link;
        }
      }

      // save col pointer of current column
      Blocks[diagBlockIndex].start_.push_back(current_nz_block);
      Blocks[linkBlockIndex].start_.push_back(current_nz_link);
    }

    assert(current_nz_block == nzCount[2 * blockId]);
    assert(current_nz_link == nzCount[2 * blockId + 1]);

    Blocks[diagBlockIndex].num_row_ = blockSize[blockId];
    Blocks[diagBlockIndex].num_col_ = blockSize[blockId];
    Blocks[linkBlockIndex].num_row_ = blockSize.back();
    Blocks[linkBlockIndex].num_col_ = blockSize[blockId];

    // print blocks for debugging
    if (debug) {
      char str[50];
      snprintf(str, 50, "debug_data/block%d.txt", diagBlockIndex);
      debug_print(Blocks[diagBlockIndex], str);
      snprintf(str, 50, "debug_data/block%d.txt", linkBlockIndex);
      debug_print(Blocks[linkBlockIndex], str);
    }

    colStart += blockSize[blockId];
  }

  // build final "Schur complement" block
  int blockIndex = 2 * nparts;

  // allocate space for block
  Blocks[blockIndex].start_.reserve(blockSize.back() + 1);
  Blocks[blockIndex].index_.reserve(nzCount[2 * nparts + 1]);
  Blocks[blockIndex].value_.reserve(nzCount[2 * nparts + 1]);

  int current_nz_schur{};
  for (int i = colStart; i < colStart + blockSize.back(); ++i) {
    int col = permutation[i];

    if (type == kMetisAugmented) {
      // diagonal is not included in augmented system
      Blocks[blockIndex].index_.push_back(i - colStart);
      Blocks[blockIndex].value_.push_back(1.0);
      ++current_nz_schur;
    }

    // go through the column
    for (int colEl = M.start_[col]; colEl < M.start_[col + 1]; ++colEl) {
      int permuted = perminv[M.index_[colEl]];

      // determine if current element belongs to diagonal block or linking
      // block
      if (permuted >= threshold) {
        Blocks[blockIndex].index_.push_back(permuted - threshold);
        Blocks[blockIndex].value_.push_back(M.value_[colEl]);
        ++current_nz_schur;
      }
    }

    // save col pointer of current column
    Blocks[blockIndex].start_.push_back(current_nz_schur);
  }

  assert(current_nz_schur == nzCount[2 * nparts]);

  Blocks[blockIndex].num_row_ = blockSize.back();
  Blocks[blockIndex].num_col_ = blockSize.back();

  // print block for debugging
  if (debug) {
    char str[50];
    snprintf(str, 50, "debug_data/block%d.txt", blockIndex);
    debug_print(Blocks[blockIndex], str);
  }
}

void Metis_caller::getNonzeros() {
  // There are nparts diagonal blocks, nparts linking blocks and 1 Schur block.
  // nzCount has 2 * (nparts + 1) entries.
  // For 0 <= i < nparts:
  // - nzCount[2 * i] is the number of nonzeros in diagonal block i
  // - nzCount[2 * i + 1] is the number of nonzeros in linking block i
  // - nzCount[2 * nparts] is the number of nonzeros in the Schur
  //    block.
  // - nzCount[2 * nparts + 1] does non represent any real block, it is
  //    used as a sum check at the end.

  assert(nzCount.size() == blockSize.size() * 2);

  // get cumulative sum of blockSize
  std::vector<int> thresholds(blockSize);
  for (int i = 1; i < thresholds.size(); ++i) {
    thresholds[i] += thresholds[i - 1];
  }

  // go through the nodes
  for (int node = 0; node < M.num_row_; ++node) {
    int partNode = partition[node];

    // go through the neighbours
    for (int j = M.start_[node]; j < M.start_[node + 1]; ++j) {
      int neigh = M.index_[j];

      // skip self loops (diagonal nonzeros)
      if (neigh == node) continue;

      // count one nonzero in the right position
      if (partNode == partition[neigh]) {
        ++nzCount[2 * partNode];
      } else {
        ++nzCount[2 * partNode + 1];
      }
    }
  }

  // add diagonal nonzeros & check
  int check{};
  for (int i = 0; i < blockSize.size(); ++i) {
    nzCount[2 * i] += blockSize[i];
    check += nzCount[2 * i + 1];
  }
  assert(check == 2 * nzCount.back());
}

void Metis_caller::debug_print(const HighsSparseMatrix& mat,
                               const std::string& filename) {
  out_file.open(filename);
  out_file << mat.start_.size() - 1 << '\n';
  for (int i : mat.start_) {
    out_file << i << '\n';
  }
  for (int i : mat.index_) {
    out_file << i << '\n';
  }
  for (double i : mat.value_) {
    out_file << i << '\n';
  }
  out_file.close();
}
