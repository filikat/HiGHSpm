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

void getMetisPermutation(const HighsSparseMatrix& A, int type, int nparts,
                         std::vector<int>& permutation) {
  idx_t nvertex{};
  idx_t nedges{};
  std::vector<idx_t> adj_ptr;
  std::vector<idx_t> adj_lst;

  // -----------------------------------------------------------
  // set up the augmented system
  // -----------------------------------------------------------
  if (type == kMetisAugmented) {
    nvertex = A.num_row_ + A.num_col_;
    nedges = A.numNz() * 2;

    // allocate space for augmented matrix
    adj_ptr.resize(nvertex + 1);
    adj_lst.resize(nedges);

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
    }
    int adj_shift = A.numNz();
    for (int i = 0; i < A.numNz(); ++i) {
      adj_lst[i + adj_shift] = A_t.index_[i];
    }
  }
  // -----------------------------------------------------------
  // set up the normal equations
  // -----------------------------------------------------------
  else if (type == kMetisNormalEq) {
    HighsSparseMatrix AAt;
    std::vector<double> theta(A.num_col_, 1.0);
    computeAThetaAT(A, theta, AAt);
    adj_ptr = AAt.start_;
    adj_lst = AAt.index_;
    nvertex = A.num_row_;
    nedges = AAt.numNz();
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

  callMetis(nvertex, 1, adj_ptr.data(), adj_lst.data(), nparts, options,
            &objval, partition.data());
  // -----------------------------------------------------------

  // -----------------------------------------------------------
  // find permutations based on partition
  // -----------------------------------------------------------
  std::vector<int> permutationMM(nvertex);
  std::vector<int> blockSizeMM(nparts + 1);
  vertexCoverMM(nvertex, nedges, nparts, partition, adj_ptr, adj_lst,
                permutationMM, blockSizeMM);

  std::vector<int> permutationN(nvertex);
  std::vector<int> blockSizeN(nparts + 1);
  vertexCoverN(nvertex, nedges, nparts, partition, adj_ptr, adj_lst,
               permutationN, blockSizeN);

  // select permutation with smallest Schur complement
  if (blockSizeMM.back() > blockSizeN.back()) {
    permutation = permutationN;
  } else {
    permutation = permutationMM;
  }
  // -----------------------------------------------------------

  // -----------------------------------------------------------
  // save to file for debugging
  // -----------------------------------------------------------
  if (false) {
    std::ofstream out_file;

    out_file.open("partition.txt");
    for (int i : partition) {
      out_file << i << '\n';
    }
    out_file.close();

    out_file.open("permMM.txt");
    for (int i : permutationMM) {
      out_file << i << '\n';
    }
    out_file.close();

    out_file.open("permN.txt");
    for (int i : permutationN) {
      out_file << i << '\n';
    }
    out_file.close();

    out_file.open("blockSizeMM.txt");
    for (int i : blockSizeMM) {
      out_file << i << '\n';
    }
    out_file.close();

    out_file.open("blockSizeN.txt");
    for (int i : blockSizeN) {
      out_file << i << '\n';
    }
    out_file.close();

    out_file.open("A_ptr.txt");
    for (int i : A.start_) {
      out_file << i << '\n';
    }
    out_file.close();

    out_file.open("A_adj.txt");
    for (int i : A.index_) {
      out_file << i << '\n';
    }
    out_file.close();
  }
  // -----------------------------------------------------------
}
