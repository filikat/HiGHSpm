#include "Metis_caller.h"

#include <fstream>

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

void getMetisPermutation(const HighsSparseMatrix& A, MetisPartitionType type,
                         int nparts, std::vector<int>& permutation,
                         std::vector<int>& blockSize) {
  idx_t nvertex{};
  idx_t nedges{};

  HighsSparseMatrix M;

  // -----------------------------------------------------------
  // set up the augmented system
  // -----------------------------------------------------------
  if (type == kMetisAugmented) {
    nvertex = A.num_row_ + A.num_col_;
    nedges = A.numNz() * 2;

    // allocate space for augmented matrix
    std::vector<idx_t> adj_ptr(nvertex + 1);
    std::vector<idx_t> adj_lst(nedges);
    std::vector<double> adj_val(nedges);

    // temporary A transpose
    HighsSparseMatrix A_t = A;
    A_t.ensureRowwise();

    // create pointers of augmented matrix
    for (int i = 0; i < A.num_col_ + 1; ++i) {
      adj_ptr[i] = A.start_[i];
    }
    int shift = A.num_col_;
    int ptr_shift = A.numNz();
    for (int i = 1; i < A.num_row_ + 1; ++i) {
      adj_ptr[i + shift] = A_t.start_[i] + ptr_shift;
    }

    // create adjacency list of augmented matrix
    for (int i = 0; i < A.numNz(); ++i) {
      adj_lst[i] = A.index_[i] + A.num_col_;
      adj_val[i] = A.value_[i];
    }
    int adj_shift = A.numNz();
    for (int i = 0; i < A.numNz(); ++i) {
      adj_lst[i + adj_shift] = A_t.index_[i];
      adj_val[i + adj_shift] = A_t.value_[i];
    }

    M.num_col_ = nvertex;
    M.num_row_ = nvertex;
    M.start_ = adj_ptr;
    M.index_ = adj_lst;
    M.value_ = adj_val;
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
  // set up and call metis
  // -----------------------------------------------------------
  // create space for outputs
  idx_t objval{};
  std::vector<idx_t> partition(nvertex);

  // initialize metis options
  idx_t options[METIS_NOPTIONS];
  setOptions(options);

  callMetis(nvertex, 1, M.start_.data(), M.index_.data(), nparts, options,
            &objval, partition.data());
  // -----------------------------------------------------------

  // -----------------------------------------------------------
  // find permutations based on partition
  // -----------------------------------------------------------
  std::vector<int> permutationMM(nvertex);
  std::vector<int> blockSizeMM(nparts + 1);
  vertexCoverMM(nvertex, nedges, nparts, partition, M.start_, M.index_,
                permutationMM, blockSizeMM);

  std::vector<int> permutationN(nvertex);
  std::vector<int> blockSizeN(nparts + 1);
  vertexCoverN(nvertex, nedges, nparts, partition, M.start_, M.index_,
               permutationN, blockSizeN);

  // select permutation with smallest Schur complement
  if (blockSizeMM.back() > blockSizeN.back()) {
    permutation = permutationN;
    blockSize = blockSizeN;
  } else {
    permutation = permutationMM;
    blockSize = blockSizeMM;
  }
  // -----------------------------------------------------------

  std::vector<HighsSparseMatrix> Blocks(2 * nparts + 1);
  getBlocks(permutation, blockSize, type, M, Blocks);

  // -----------------------------------------------------------
  // save to file for debugging
  // -----------------------------------------------------------
  if (true) {
    std::ofstream out_file;

    out_file.open("debug_data/partition.txt");
    for (int i : partition) {
      out_file << i << '\n';
    }
    out_file.close();

    out_file.open("debug_data/permMM.txt");
    for (int i : permutationMM) {
      out_file << i << '\n';
    }
    out_file.close();

    out_file.open("debug_data/permN.txt");
    for (int i : permutationN) {
      out_file << i << '\n';
    }
    out_file.close();

    out_file.open("debug_data/blockSizeMM.txt");
    for (int i : blockSizeMM) {
      out_file << i << '\n';
    }
    out_file.close();

    out_file.open("debug_data/blockSizeN.txt");
    for (int i : blockSizeN) {
      out_file << i << '\n';
    }
    out_file.close();

    out_file.open("debug_data/A_ptr.txt");
    for (int i : A.start_) {
      out_file << i << '\n';
    }
    out_file.close();

    out_file.open("debug_data/A_adj.txt");
    for (int i : A.index_) {
      out_file << i << '\n';
    }
    out_file.close();

    out_file.open("debug_data/A_val.txt");
    for (double i : A.value_) {
      out_file << i << '\n';
    }
    out_file.close();
  }
  // -----------------------------------------------------------
}

void getBlocks(const std::vector<int>& permutation,
               const std::vector<int>& blockSize, MetisPartitionType type,
               const HighsSparseMatrix& M,
               std::vector<HighsSparseMatrix>& Blocks) {
  // get inverse permutation
  std::vector<int> perminv(permutation.size());
  for (int i = 0; i < perminv.size(); ++i) {
    perminv[permutation[i]] = i;
  }

  int threshold = M.num_col_ - blockSize.back();
  int nparts = blockSize.size() - 1;

  int colStart = 0;

  std::ofstream out_file;

  // go through the blocks
  for (int blockId = 0; blockId < nparts; ++blockId) {
    // indices to access the correct block in Blocks
    int diagBlockIndex = 2 * blockId;
    int linkBlockIndex = 2 * blockId + 1;

    // count nonzeros in diagonal block and linking block
    int current_nz_block{};
    int current_nz_link{};

    Blocks[diagBlockIndex].start_.resize(blockSize[blockId] + 1);
    Blocks[linkBlockIndex].start_.resize(blockSize[blockId] + 1);

    // TO DO preallocate index_ and value_  of Blocks in some way

    // go through the columns in the order of permutation
    for (int i = colStart; i < colStart + blockSize[blockId]; ++i) {
      int col = permutation[i];

      Blocks[diagBlockIndex].start_[i - colStart] = (current_nz_block);
      Blocks[linkBlockIndex].start_[i - colStart] = (current_nz_link);

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
    }
    Blocks[diagBlockIndex].start_.back() = (current_nz_block);
    Blocks[linkBlockIndex].start_.back() = (current_nz_link);

    Blocks[diagBlockIndex].num_row_ = blockSize[blockId];
    Blocks[diagBlockIndex].num_col_ = blockSize[blockId];
    Blocks[linkBlockIndex].num_row_ = blockSize.back();
    Blocks[linkBlockIndex].num_col_ = blockSize[blockId];

    // -----------------------------------------------------------
    // Print blocks for debugging
    // -----------------------------------------------------------
    char str[50];

    snprintf(str, 50, "debug_data/block%d.txt", diagBlockIndex);
    out_file.open(str);
    out_file << Blocks[diagBlockIndex].start_.size() - 1 << '\n';
    for (int i : Blocks[diagBlockIndex].start_) {
      out_file << i << '\n';
    }
    for (int i : Blocks[diagBlockIndex].index_) {
      out_file << i << '\n';
    }
    for (double i : Blocks[diagBlockIndex].value_) {
      out_file << i << '\n';
    }
    out_file.close();

    snprintf(str, 50, "debug_data/block%d.txt", linkBlockIndex);
    out_file.open(str);
    out_file << Blocks[linkBlockIndex].start_.size() - 1 << '\n';
    for (int i : Blocks[linkBlockIndex].start_) {
      out_file << i << '\n';
    }
    for (int i : Blocks[linkBlockIndex].index_) {
      out_file << i << '\n';
    }
    for (double i : Blocks[linkBlockIndex].value_) {
      out_file << i << '\n';
    }
    out_file.close();
    // -----------------------------------------------------------

    colStart += blockSize[blockId];
  }

  // build final "Schur complement" block
  int blockIndex = 2 * nparts;
  Blocks[blockIndex].start_.resize(blockSize.back() + 1);

  int current_nz_schur{};
  for (int i = colStart; i < colStart + blockSize.back(); ++i) {
    int col = permutation[i];
    Blocks[blockIndex].start_[i - colStart] = (current_nz_schur);
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
  }
  Blocks[blockIndex].start_.back() = (current_nz_schur);

  Blocks[blockIndex].num_row_ = blockSize.back();
  Blocks[blockIndex].num_col_ = blockSize.back();

  // -----------------------------------------------------------
  // Print block for debugging
  // -----------------------------------------------------------
  char str[50];

  snprintf(str, 50, "debug_data/block%d.txt", blockIndex);
  out_file.open(str);
  out_file << Blocks[blockIndex].start_.size() - 1 << '\n';
  for (int i : Blocks[blockIndex].start_) {
    out_file << i << '\n';
  }
  for (int i : Blocks[blockIndex].index_) {
    out_file << i << '\n';
  }
  for (double i : Blocks[blockIndex].value_) {
    out_file << i << '\n';
  }
  out_file.close();
  // -----------------------------------------------------------
}