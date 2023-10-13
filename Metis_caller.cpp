#include "Metis_caller.h"

// -----------------------------------------------------------
// Metis wrapper
// -----------------------------------------------------------
// extern "C" {
void metis_wrapper_set_options(idx_t* options) {
  METIS_SetDefaultOptions(options);
}
void metis_wrapper_call_metis(idx_t nvertex, idx_t nconstraints, idx_t* adj_ptr,
                              idx_t* adj_lst, idx_t nparts, idx_t* options,
                              idx_t* objval, idx_t* part) {
  idx_t status =
      METIS_PartGraphKway(&nvertex, &nconstraints, adj_ptr, adj_lst, NULL, NULL,
                          NULL, &nparts, NULL, NULL, options, objval, part);

  assert(status == METIS_OK);
}
//}
// -----------------------------------------------------------

Metis_caller::Metis_caller(const HighsSparseMatrix& input_A,
                           MetisPartitionType input_type, int input_nparts) {
  nparts = input_nparts;
  type = input_type;
  A = &input_A;

  // -----------------------------------------------------------
  // set up the augmented system
  // -----------------------------------------------------------
  if (type == kMetisAugmented) {
    nvertex = A->num_row_ + A->num_col_;
    nedges = A->numNz() * 2;

    // allocate space for augmented matrix
    M.num_col_ = nvertex;
    M.num_row_ = nvertex;
    M.start_.resize(nvertex + 1);
    M.index_.resize(nedges);
    M.value_.resize(nedges);

    // temporary A transpose
    HighsSparseMatrix A_t = *A;
    A_t.ensureRowwise();

    // create pointers of augmented matrix
    for (int i = 0; i < A->num_col_ + 1; ++i) {
      M.start_[i] = A->start_[i];
    }
    int shift = A->num_col_;
    int ptr_shift = A->numNz();
    for (int i = 1; i < A->num_row_ + 1; ++i) {
      M.start_[i + shift] = A_t.start_[i] + ptr_shift;
    }

    // create adjacency list of augmented matrix
    for (int i = 0; i < A->numNz(); ++i) {
      M.index_[i] = A->index_[i] + A->num_col_;
      M.value_[i] = A->value_[i];
    }
    int adj_shift = A->numNz();
    for (int i = 0; i < A->numNz(); ++i) {
      M.index_[i + adj_shift] = A_t.index_[i];
      M.value_[i + adj_shift] = A_t.value_[i];
    }
  }
  // -----------------------------------------------------------
  // set up the normal equations
  // -----------------------------------------------------------
  else if (type == kMetisNormalEq) {
    std::vector<double> theta(A->num_col_, 1.0);
    computeAThetaAT(*A, theta, M);
    nvertex = A->num_row_;
    nedges = M.numNz();
  } else {
    std::cerr << "Wront type of matrix for Metis parition\n";
    return;
  }
}

void Metis_caller::getPartition() {
  // create space for outputs
  idx_t objval{};
  partition.resize(nvertex);

  // initialize metis options
  idx_t options[METIS_NOPTIONS];
  metis_wrapper_set_options(options);

  // call Metis to get the partition
  metis_wrapper_call_metis(nvertex, 1, M.start_.data(), M.index_.data(), nparts,
                           options, &objval, partition.data());
}

void Metis_caller::getPermutation() {
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
    debug_print(A->start_, "debug_data/A_ptr.txt");
    debug_print(A->index_, "debug_data/A_adj.txt");
    debug_print(A->value_, "debug_data/A_val.txt");
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

  // get inverse permutation
  perminv.resize(nvertex);
  for (int i = 0; i < perminv.size(); ++i) {
    perminv[permutation[i]] = i;
  }
}

void Metis_caller::getBlocks(const std::vector<double>& diag1,
                             const std::vector<double>& diag2) {
  // normal equations has to be recomputed with correct diagonal
  // (cannot be easily updated)
  if (type == kMetisNormalEq) {
    M.clear();
    computeAThetaAT(*A, diag1, M);
  }

  // if Blocks were already computed, aug system can be easily updated
  if (!Blocks.empty() && type == kMetisAugmented) {
    updateDiag(diag1, diag2);
    return;
  }

  int threshold = M.num_col_ - blockSize.back();

  // get number of nonzeros in blocks for preallocation
  // (just the first time)
  if (nzCount.empty()) getNonzeros();

  // allocate/clear space for blocks
  Blocks.assign(2 * nparts + 1, HighsSparseMatrix());

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

        // extract diagonal element from diag1 or diag2
        double diagEl =
            col < diag1.size() ? diag1[col] : diag2[col - diag1.size()];
        Blocks[diagBlockIndex].value_.push_back(diagEl);
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
  Blocks[blockIndex].index_.reserve(nzCount[2 * nparts]);
  Blocks[blockIndex].value_.reserve(nzCount[2 * nparts]);

  int current_nz_schur{};
  for (int i = colStart; i < colStart + blockSize.back(); ++i) {
    int col = permutation[i];

    if (type == kMetisAugmented) {
      // diagonal is not included in augmented system
      Blocks[blockIndex].index_.push_back(i - colStart);
      // extract diagonal element from diag1 or diag2
      double diagEl =
          col < diag1.size() ? diag1[col] : diag2[col - diag1.size()];
      Blocks[blockIndex].value_.push_back(diagEl);
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

  nzCount.resize(2 * nparts + 2);

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

void Metis_caller::updateDiag(const std::vector<double>& diag1,
                              const std::vector<double>& diag2) {
  assert(type == kMetisAugmented);

  // index to access permutation
  int permIndex{};

  // go through the diagonal blocks
  for (int partId = 0; partId <= nparts; ++partId) {
    HighsSparseMatrix* curBlock = &Blocks[2 * partId];

    // go through the columns
    for (int col = 0; col < curBlock->num_col_; ++col) {
      int firstEl = curBlock->start_[col];

      // first element of the column is diagonal element
      assert(curBlock->index_[firstEl] == col);

      // which element of the original diagonal corresponds to the current
      // column
      int diagIndex = permutation[permIndex];

      // extract the correct diagonal element
      double newDiagEl = diagIndex < diag1.size()
                             ? diag1[diagIndex]
                             : diag2[diagIndex - diag1.size()];

      curBlock->value_[firstEl] = newDiagEl;

      ++permIndex;
    }

    if (debug) {
      char str[50];
      snprintf(str, 50, "debug_data/block%d.txt", 2 * partId);
      debug_print(*curBlock, str);
    }
  }
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

void Metis_caller::factor() {
  // allocate space for factorizations of blocks
  invertData.resize(nparts + 1);
  expData.resize(nparts + 1);

  // to build Schur complement:
  // for each linking block
  // - access rows of linking block (stored row-wise)
  // - create dense vector with row
  // - solve with diagonal block factorization
  // - product with linking block (stored col-wise)
  // - obtain column of Schur complement
  // Store everything in a dense matrix

  auto t1 = getWallTime();

  int linkSize = blockSize.back();

  // space for dense schur complement
  std::vector<std::vector<double>> schurComplement(
      linkSize, std::vector<double>(linkSize, 0.0));

  // insert original schur block
  HighsSparseMatrix& schurB = Blocks[2 * nparts];
  for (int col = 0; col < linkSize; ++col) {
    for (int el = schurB.start_[col]; el < schurB.start_[col + 1]; ++el) {
      schurComplement[col][schurB.index_[el]] += schurB.value_[el];
    }
  }

  // go through the parts
  for (int i = 0; i < nparts; ++i) {
    // factorize the diagonal blocks
    expData[i].reset();
    blockInvert(Blocks[2 * i], invertData[i], expData[i]);

    // current linking block
    HighsSparseMatrix& B = Blocks[2 * i + 1];

    // row-wise linking block
    HighsSparseMatrix Bt = B;
    Bt.ensureRowwise();

    // go through the rows of B
    for (int row = 0; row < B.num_row_; ++row) {
      // space for dense vector
      std::vector<double> denseRow(B.num_col_, 0.0);
      // space for column of the schur complement
      std::vector<double> schurContribution(linkSize, 0.0);

      // go through the entries of the row and fill in dense_row
      for (int elem = Bt.start_[row]; elem < Bt.start_[row + 1]; ++elem) {
        denseRow[Bt.index_[elem]] = Bt.value_[elem];
      }

      // solve with diagonal block
      blockSolve(denseRow.data(), 1, invertData[i], expData[i]);

      // multiply by B
      B.product(schurContribution, denseRow);

      // subtract contribution of current row from schur complement
      denseMatrixPlusVector(schurComplement, schurContribution, row, -1.0);
    }
  }
  std::cout << "time " << getWallTime() - t1 << '\n';

  // convert schur complement to sparse matrix
  HighsSparseMatrix sparseSchur;
  sparseSchur.start_.reserve(linkSize + 1);
  sparseSchur.index_.reserve(linkSize * linkSize);
  sparseSchur.value_.reserve(linkSize * linkSize);
  for (int col = 0; col < linkSize; ++col) {
    sparseSchur.start_.push_back(sparseSchur.start_[col] + linkSize);
    for (int row = 0; row < linkSize; ++row) {
      sparseSchur.index_.push_back(row);
      sparseSchur.value_.push_back(schurComplement[col][row]);
    }
  }
  sparseSchur.num_row_ = linkSize;
  sparseSchur.num_col_ = linkSize;
  schurComplement.clear();

  // factorize schur complement
  blockInvert(sparseSchur, invertData.back(), expData.back());
}

void Metis_caller::solve(const std::vector<double>& rhs,
                         std::vector<double>& lhs) {
  // space for permuted rhs and lhs
  std::vector<double> perm_rhs(rhs.size());
  std::vector<double> perm_lhs(lhs.size());

  // permute rhs
  for (int i = 0; i < rhs.size(); ++i) {
    perm_rhs[i] = rhs[permutation[i]];
  }

  // allocate space for partitioning rhs and lhs
  block_rhs.resize(nparts + 1);
  block_lhs.resize(nparts + 1);

  // rhs of schur complement is temporarily stored in the last solution block
  block_lhs.back().resize(blockSize.back());
  std::vector<double>& schur_rhs(block_lhs.back());

  int start{};
  for (int i = 0; i <= nparts; ++i) {
    // allocate space for current rhs and lhs blocks
    block_rhs[i].resize(blockSize[i]);
    block_lhs[i].resize(blockSize[i]);

    // break rhs into blocks
    for (int j = 0; j < blockSize[i]; ++j) {
      block_rhs[i][j] = perm_rhs[start + j];
    }
    start += blockSize[i];

    // contribution to schur_rhs:
    // schur_rhs = rhs_link - \sum_i Bi * Di^-1 * rhs_i
    if (i < nparts) {
      std::vector<double> temp(block_rhs[i]);
      std::vector<double> schur_rhs_temp(blockSize.back());
      blockSolve(temp.data(), 1, invertData[i], expData[i]);
      Blocks[2 * i + 1].product(schur_rhs_temp, temp);
      VectorAdd(schur_rhs, schur_rhs_temp, -1.0);
    } else {
      VectorAdd(schur_rhs, block_rhs[i], 1.0);
    }
  }

  // solve linear system with schur complement to compute last block of solution
  blockSolve(schur_rhs.data(), 1, invertData.back(), expData.back());

  // compute the other blocks of the solution:
  // lhs_i = Di^-1 * (rhs_i - Bi^T * lhs_last)
  std::vector<double>& lastSol(block_lhs.back());
  int positionInLhs{};
  for (int i = 0; i < nparts; ++i) {
    Blocks[2 * i + 1].productTranspose(block_lhs[i], lastSol);
    VectorAdd(block_lhs[i], block_rhs[i], -1.0);
    VectorScale(block_lhs[i], -1.0);
    blockSolve(block_lhs[i].data(), 1, invertData[i], expData[i]);

    // copy solution into lhs
    std::copy(block_lhs[i].begin(), block_lhs[i].end(),
              perm_lhs.begin() + positionInLhs);
    positionInLhs += blockSize[i];
  }

  // copy last block of solution
  std::copy(lastSol.begin(), lastSol.end(), perm_lhs.begin() + positionInLhs);

  // permute back lhs
  for (int i = 0; i < lhs.size(); ++i) {
    lhs[i] = perm_lhs[perminv[i]];
  }

  debug_print(rhs, "debug_data/rhs.txt");
  debug_print(lhs, "debug_data/lhs.txt");
}

void denseMatrixPlusVector(std::vector<std::vector<double>>& mat,
                           const std::vector<double>& vec, int col,
                           double alpha) {
  // update column col of dense matrix mat by adding vector vec scaled by alpha
  for (int i = 0; i < vec.size(); ++i) {
    mat[col][i] += alpha * vec[i];
  }
}
