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

Metis_caller::Metis_caller(const HighsSparseMatrix& input_A, int input_type,
                           int input_nparts) {
  nparts = input_nparts;
  type = input_type;
  A = &input_A;

  double t0 = getWallTime();

  // -----------------------------------------------------------
  // set up the augmented system
  // -----------------------------------------------------------
  if (type == kOptionNlaMetisAugmented) {
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
  else if (type == kOptionNlaMetisNormalEq) {
    std::vector<double> theta(A->num_col_, 1.0);
    // To use metis, what matters is the sparsity pattern of M. Compute M using
    // temporary A, with all nonzeros equal to 1, to avoid numerical
    // cancellation.
    HighsSparseMatrix tempA = *A;
    tempA.value_.assign(A->value_.size(), 1.0);
    computeAThetaAT(tempA, theta, M);
    nvertex = A->num_row_;
    nedges = M.numNz();
  } else {
    std::cerr << "Wrong type of matrix for Metis partition\n";
    return;
  }

  // space for factorizations
  invertData.resize(nparts + 1, {});
  expData.resize(nparts + 1);

  initial_time += getWallTime() - t0;
}

void Metis_caller::getPartition() {
  double t0 = getWallTime();

  // create space for outputs
  idx_t objval{};
  partition.resize(nvertex);

  // initialize metis options
  idx_t options[METIS_NOPTIONS];
  metis_wrapper_set_options(options);

  // call Metis to get the partition
  metis_wrapper_call_metis(nvertex, 1, M.start_.data(), M.index_.data(), nparts,
                           options, &objval, partition.data());

  initial_time += getWallTime() - t0;
}

void Metis_caller::getPermutation() {
  double t0 = getWallTime();

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
  if (debug) {
    debug_print(partition, "debug_data/partition.txt");
  }

  // get inverse permutation
  perminv.resize(nvertex);
  for (int i = 0; i < perminv.size(); ++i) {
    perminv[permutation[i]] = i;
  }

  initial_time += getWallTime() - t0;

  // get number of nonzeros in blocks for preallocation
  getNonzeros();
}

void Metis_caller::getBlocks(const std::vector<double>& diag1,
                             const std::vector<double>& diag2) {
  double t0 = getWallTime();

  // normal equations has to be recomputed with correct diagonal
  // (cannot be easily updated)
  if (type == kOptionNlaMetisNormalEq) {
    M.clear();
    computeAThetaAT(*A, diag1, M);
  }

  // if Blocks were already computed, aug system can be easily updated
  if (!Blocks.empty() && type == kOptionNlaMetisAugmented) {
    updateDiag(diag1, diag2);
    getBlocks_time += getWallTime() - t0;
    return;
  }

  int threshold = M.num_col_ - blockSize.back();

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

      if (type == kOptionNlaMetisAugmented) {
        // diagonal is not included in augmented system
        Blocks[diagBlockIndex].index_.push_back(i - colStart);

        // extract diagonal element from diag1 or diag2
        double diagEl =
            col < diag1.size() ? -diag1[col] : diag2[col - diag1.size()];
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

    // These assertions with equality may fail because elements of the normal
    // equations that are too small (<1e-10) are ignored, changing the actual
    // number of nonzeros depending on the values of theta.
    assert(current_nz_block <= nzCount[2 * blockId]);
    assert(current_nz_link <= nzCount[2 * blockId + 1]);

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

    if (type == kOptionNlaMetisAugmented) {
      // diagonal is not included in augmented system
      Blocks[blockIndex].index_.push_back(i - colStart);
      // extract diagonal element from diag1 or diag2
      double diagEl =
          col < diag1.size() ? -diag1[col] : diag2[col - diag1.size()];
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

  assert(current_nz_schur <= nzCount[2 * nparts]);

  Blocks[blockIndex].num_row_ = blockSize.back();
  Blocks[blockIndex].num_col_ = blockSize.back();

  // print block for debugging
  if (debug) {
    char str[50];
    snprintf(str, 50, "debug_data/block%d.txt", blockIndex);
    debug_print(Blocks[blockIndex], str);
  }

  if (debug) {
    debug_print(diag1, "debug_data/diag1.txt");
    debug_print(diag2, "debug_data/diag2.txt");
  }

  getBlocks_time += getWallTime() - t0;
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
  assert(type == kOptionNlaMetisAugmented);

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
                             ? -diag1[diagIndex]
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
  // to build Schur complement:
  // for each linking block
  // - access rows of linking block B^T (stored row-wise)
  // - create dense vector with row b
  // - do a forward solve and store as sparse column
  //    L^-1 * b
  // - do a diagonal solve and store as sparse column
  //    D^-1 * L^-1 * b
  // - leave L^-1 * B^T col-wise, but get D^-1 * L^-1 * B^T row-wise
  // - use same trick as computeAThetaAT to find contribution to Schur
  //    complement
  // - Store everything in a dense matrix

  double t0 = getWallTime();

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

  schur_initial_time += getWallTime() - t0;

  // go through the parts
  for (int i = 0; i < nparts; ++i) {
    // factorize the diagonal blocks
    double t1 = getWallTime();
    expData[i].reset();
    blockInvert(Blocks[2 * i], invertData[i], expData[i]);
    factorBlocks_time += getWallTime() - t1;
    schur_factor_time += getWallTime() - t1;

    // current linking block
    HighsSparseMatrix& B = Blocks[2 * i + 1];

    // row-wise linking block
    HighsSparseMatrix Bt = B;
    Bt.ensureRowwise();

    // space for L^-1 B^T and D^-1 L^-1 B^T
    HighsSparseMatrix forwardSolvedBlock;
    HighsSparseMatrix diagForwardSolvedBlock;
    forwardSolvedBlock.num_row_ = B.num_col_;
    forwardSolvedBlock.num_col_ = B.num_row_;
    diagForwardSolvedBlock.num_row_ = B.num_col_;
    diagForwardSolvedBlock.num_col_ = B.num_row_;

    // go through the rows of B
    for (int row = 0; row < B.num_row_; ++row) {
      t1 = getWallTime();
      // space for dense vector
      std::vector<double> denseRow(B.num_col_, 0.0);

      // avoid computations if row is empty
      if (Bt.start_[row] == Bt.start_[row + 1]) {
        forwardSolvedBlock.start_.push_back(forwardSolvedBlock.start_.back());
        diagForwardSolvedBlock.start_.push_back(
            diagForwardSolvedBlock.start_.back());
        continue;
      }

      // go through the entries of the row and fill in dense_row
      for (int elem = Bt.start_[row]; elem < Bt.start_[row + 1]; ++elem) {
        denseRow[Bt.index_[elem]] = Bt.value_[elem];
      }
      schur_fillrow_time += getWallTime() - t1;

      t1 = getWallTime();
      // on return from diagonalForwardSolve, denseRow contains L^-1 * denseRow
      // and diagForwardSolvedRow contains D^-1 * L^-1 * denseRow
      std::vector<double> diagForwardSolvedRow(B.num_col_, 0.0);
      diagonalForwardSolve(denseRow.data(), invertData[i], expData[i],
                           diagForwardSolvedRow.data());
      schur_dfsolve_time += getWallTime() - t1;

      t1 = getWallTime();
      // transform denseRow into sparse column of forwardSolvedBlock
      int countNnzRow{};
      for (int i = 0; i < denseRow.size(); ++i) {
        if (std::fabs(denseRow[i]) > 1e-20) {
          forwardSolvedBlock.index_.push_back(i);
          forwardSolvedBlock.value_.push_back(denseRow[i]);
          ++countNnzRow;
        }
      }
      forwardSolvedBlock.start_.push_back(forwardSolvedBlock.start_.back() +
                                          countNnzRow);

      // transform diagForwardSolvedRow into sparse column of
      // diagForwardSolvedBlock
      countNnzRow = 0;
      for (int i = 0; i < diagForwardSolvedRow.size(); ++i) {
        if (std::fabs(diagForwardSolvedRow[i]) > 1e-20) {
          diagForwardSolvedBlock.index_.push_back(i);
          diagForwardSolvedBlock.value_.push_back(diagForwardSolvedRow[i]);
          ++countNnzRow;
        }
      }
      diagForwardSolvedBlock.start_.push_back(
          diagForwardSolvedBlock.start_.back() + countNnzRow);

      schur_transform_time += getWallTime() - t1;
    }

    // diagForwardSolvedBlock needs row-wise access
    diagForwardSolvedBlock.ensureRowwise();

    // similar trick as computeAThetaAT, with simplifications, to compute
    //      Q^T        *          R
    // (L^-1 * B^T)^T  *  (D^-1 * L^-1 * B^T)
    // ^^^^^^^^^^^^       ^^^^^^^^^^^^^^^^^^^
    //  col-wise               row-wise
    // access to Q            access to R
    //
    // Q and R may have different sparsity pattern, because D may have 2x2
    // pivots. This may cause inefficiencies.

    HighsSparseMatrix& Q = forwardSolvedBlock;
    HighsSparseMatrix& R = diagForwardSolvedBlock;

    t1 = getWallTime();

    // go through columns of Q
    for (int col = 0; col < Q.num_col_; ++col) {
      // go through entries of the column
      for (int colEl = Q.start_[col]; colEl < Q.start_[col + 1]; ++colEl) {
        // row index of current entry
        int row = Q.index_[colEl];
        // go through entries of the row of R
        for (int rowEl = R.start_[row]; rowEl < R.start_[row + 1]; ++rowEl) {
          // this entry contributes to Schur(col,col1) and Schur (col1,col)
          int col1 = R.index_[rowEl];

          // avoid counting entries twice
          if (col1 < col) continue;

          double value = Q.value_[colEl] * R.value_[rowEl];
          schurComplement[col][col1] -= value;
          if (col != col1) {
            schurComplement[col1][col] -= value;
          }
        }
      }
    }
    schur_multiply_time += getWallTime() - t1;
  }

  formSchur_time += getWallTime() - t0;

  t0 = getWallTime();

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

  factorSchur_time += getWallTime() - t0;
}

void Metis_caller::solve(const std::vector<double>& rhs,
                         std::vector<double>& lhs) {
  double t0 = getWallTime();

  // space for blocks of rhs and lhs
  std::vector<std::vector<double>> block_rhs;
  std::vector<std::vector<double>> block_lhs;

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

  if (debug) {
    debug_print(rhs, "debug_data/rhs.txt");
    debug_print(lhs, "debug_data/lhs.txt");
  }

  solve_time += getWallTime() - t0;
}

void Metis_caller::printInfo() const {
  printf("* * * * * * * * * * *\n");
  printf("Using Metis for %s\n", type == kOptionNlaMetisAugmented
                                     ? "augmented system"
                                     : "normal equations");
  printf("Nodes: %8d\nEdges: %8d\n", nvertex, nedges);
  printf("Partitioning in %d parts\n", nparts);
  printf("Diagonal blocks size: ");
  for (int i = 0; i < nparts; ++i) printf("%2d ", blockSize[i]);
  printf("\nSchur complement size: %5d (%.2f%%)\n", blockSize.back(),
         100 * (double)blockSize.back() / nvertex);
  printf("* * * * * * * * * * *\n");
}

void Metis_caller::printTimes() const {
  double sum_time = initial_time + getBlocks_time + formSchur_time +
                    factorSchur_time + solve_time;
  sum_time = sum_time > 0 ? sum_time : 1;

  printf("Metis time profile\n");
  printf("initial       %5.2f (%5.1f%% sum)\n", initial_time,
         100 * initial_time / sum_time);
  printf("get blocks    %5.2f (%5.1f%% sum)\n", getBlocks_time,
         100 * getBlocks_time / sum_time);
  printf("form Schur    %5.2f (%5.1f%% sum), including factor blocks\n",
         formSchur_time, 100 * formSchur_time / sum_time);
  printf("factor blocks %5.2f (%5.1f%% sum)\n", factorBlocks_time,
         100 * factorBlocks_time / sum_time);
  printf("factor Schur  %5.2f (%5.1f%% sum)\n", factorSchur_time,
         100 * factorSchur_time / sum_time);
  printf("solve         %5.2f (%5.1f%% sum)\n", solve_time,
         100 * solve_time / sum_time);

  sum_time = schur_initial_time + schur_factor_time + schur_fillrow_time +
             schur_dfsolve_time + schur_transform_time + schur_multiply_time;
  printf("\nSchur form profile\n");
  printf("initial       %5.2f (%5.1f%% sum)\n", schur_initial_time,
         100 * schur_initial_time / sum_time);
  printf("factor        %5.2f (%5.1f%% sum)\n", schur_factor_time,
         100 * schur_factor_time / sum_time);
  printf("fillrow       %5.2f (%5.1f%% sum)\n", schur_fillrow_time,
         100 * schur_fillrow_time / sum_time);
  printf("dfsolve       %5.2f (%5.1f%% sum)\n", schur_dfsolve_time,
         100 * schur_dfsolve_time / sum_time);
  printf("transform     %5.2f (%5.1f%% sum)\n", schur_transform_time,
         100 * schur_transform_time / sum_time);
  printf("multiply      %5.2f (%5.1f%% sum)\n", schur_multiply_time,
         100 * schur_multiply_time / sum_time);
}
