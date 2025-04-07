#include "Auxiliary.h"

#include <stack>

#include "CallAndTimeBlas.h"
#include "DataCollector.h"

void counts2Ptr(std::vector<int>& ptr, std::vector<int>& w) {
  // Given the column counts in the vector w (of size n),
  // compute the column pointers in the vector ptr (of size n+1),
  // and copy the first n pointers back into w.

  int temp_nz{};
  int n = w.size();
  for (int j = 0; j < n; ++j) {
    ptr[j] = temp_nz;
    temp_nz += w[j];
    w[j] = ptr[j];
  }
  ptr[n] = temp_nz;
}

void inversePerm(const std::vector<int>& perm, std::vector<int>& iperm) {
  // Given the permutation perm, produce the inverse permutation iperm.
  // perm[i] : i-th entry to use in the new order.
  // iperm[i]: where entry i is located in the new order.

  for (int i = 0; i < perm.size(); ++i) {
    iperm[perm[i]] = i;
  }
}

void subtreeSize(const std::vector<int>& parent, std::vector<int>& sizes) {
  // Compute sizes of subtrees of the tree given by parent

  int n = parent.size();
  sizes.assign(n, 1);

  for (int i = 0; i < n; ++i) {
    int k = parent[i];
    if (k != -1) sizes[k] += sizes[i];
  }
}

void transpose(const std::vector<int>& ptr, const std::vector<int>& rows,
               std::vector<int>& ptrT, std::vector<int>& rowsT) {
  // Compute the transpose of the matrix and return it in rowsT and ptrT

  int n = ptr.size() - 1;

  std::vector<int> work(n);

  // count the entries in each row into work
  for (int i = 0; i < ptr.back(); ++i) {
    ++work[rows[i]];
  }

  // sum row sums to obtain pointers
  counts2Ptr(ptrT, work);

  for (int j = 0; j < n; ++j) {
    for (int el = ptr[j]; el < ptr[j + 1]; ++el) {
      int i = rows[el];

      // entry (i,j) becomes entry (j,i)
      int pos = work[i]++;
      rowsT[pos] = j;
    }
  }
}

void transpose(const std::vector<int>& ptr, const std::vector<int>& rows,
               const std::vector<double>& val, std::vector<int>& ptrT,
               std::vector<int>& rowsT, std::vector<double>& valT) {
  // Compute the transpose of the matrix and return it in rowsT, ptrT and valT

  int n = ptr.size() - 1;

  std::vector<int> work(n);

  // count the entries in each row into work
  for (int i = 0; i < ptr.back(); ++i) {
    ++work[rows[i]];
  }

  // sum row sums to obtain pointers
  counts2Ptr(ptrT, work);

  for (int j = 0; j < n; ++j) {
    for (int el = ptr[j]; el < ptr[j + 1]; ++el) {
      int i = rows[el];

      // entry (i,j) becomes entry (j,i)
      int pos = work[i]++;
      rowsT[pos] = j;
      valT[pos] = val[el];
    }
  }
}

void symProduct(const std::vector<int>& ptr, const std::vector<int>& rows,
                const std::vector<double>& vals, const std::vector<double>& x,
                std::vector<double>& y, double alpha) {
  // Matrix-vector product in CSC format, for symmetric matrix which stores only
  // the lower triangle.
  // Compute y = y + alpha * M * x

  const int n = ptr.size() - 1;

  for (int col = 0; col < n; ++col) {
    for (int el = ptr[col]; el < ptr[col + 1]; ++el) {
      int row = rows[el];
      double val = vals[el];

      y[row] += alpha * val * x[col];
      if (row != col) y[col] += alpha * val * x[row];
    }
  }
}

void symProductQuad(const std::vector<int>& ptr, const std::vector<int>& rows,
                    const std::vector<double>& vals,
                    const std::vector<double>& x, std::vector<HighsCDouble>& y,
                    double alpha) {
  // Matrix-vector product in CSC format, for symmetric matrix which stores only
  // the lower triangle.
  // Compute y = y + alpha * M * x

  const int n = ptr.size() - 1;

  for (int col = 0; col < n; ++col) {
    for (int el = ptr[col]; el < ptr[col + 1]; ++el) {
      int row = rows[el];
      HighsCDouble val = vals[el];

      y[row] += val * (HighsCDouble)x[col] * (HighsCDouble)alpha;
      if (row != col)
        y[col] += val * (HighsCDouble)x[row] * (HighsCDouble)alpha;
    }
  }
}

void childrenLinkedList(const std::vector<int>& parent, std::vector<int>& head,
                        std::vector<int>& next) {
  // Create linked lists of children in elimination tree.
  // parent gives the dependencies of the tree,
  // head[node] is the first child of node,
  // next[head[node]] is the second child,
  // next[next[head[node]]] is the third child...
  // until -1 is reached.

  int n = parent.size();
  head.assign(n, -1);
  next.assign(n, -1);
  for (int node = n - 1; node >= 0; --node) {
    if (parent[node] == -1) continue;
    next[node] = head[parent[node]];
    head[parent[node]] = node;
  }
}

void reverseLinkedList(std::vector<int>& head, std::vector<int>& next) {
  // Reverse the linked list of children of each node.
  // If a node has children (a -> b -> c -> -1), the reverse list contains
  // children (c -> b -> a -> -1).

  const int n = head.size();

  for (int node = 0; node < n; ++node) {
    int prev_node = -1;
    int curr_node = head[node];
    int next_node = -1;

    while (curr_node != -1) {
      next_node = next[curr_node];
      next[curr_node] = prev_node;
      prev_node = curr_node;
      curr_node = next_node;
    }

    head[node] = prev_node;
  }
}

void dfsPostorder(int node, int& start, std::vector<int>& head,
                  const std::vector<int>& next, std::vector<int>& order) {
  // Perform depth first search starting from root node and order the nodes
  // starting from the value start. head and next contain the linked list of
  // children.

  std::stack<int> stack;
  stack.push(node);

  while (!stack.empty()) {
    const int current = stack.top();
    const int child = head[current];

    if (child == -1) {
      // no children left to order,
      // remove from the stack and order
      stack.pop();
      order[start++] = current;
    } else {
      // at least one child left to order,
      // add it to the stack and remove it from the list of children
      stack.push(child);
      head[current] = next[child];
    }
  }
}

void processEdge(int j, int i, const std::vector<int>& first,
                 std::vector<int>& maxfirst, std::vector<int>& delta,
                 std::vector<int>& prevleaf, std::vector<int>& ancestor) {
  // Process edge of skeleton matrix.
  // Taken from Tim Davis "Direct Methods for Sparse Linear Systems".

  // j not a leaf of ith row subtree
  if (i <= j || first[j] <= maxfirst[i]) {
    return;
  }

  // max first[j] so far
  maxfirst[i] = first[j];

  // previous leaf of ith row subtree
  int jprev = prevleaf[i];

  // A(i,j) is in the skeleton matrix
  delta[j]++;

  if (jprev != -1) {
    // find least common ancestor of jprev and j
    int q = jprev;
    while (q != ancestor[q]) {
      q = ancestor[q];
    }

    // path compression
    int sparent;
    for (int s = jprev; s != q; s = sparent) {
      sparent = ancestor[s];
      ancestor[s] = q;
    }

    // consider overlap
    delta[q]--;
  }

  // previous leaf of ith subtree set to j
  prevleaf[i] = j;
}

double getDiagStart(int n, int k, int nb, int n_blocks, std::vector<int>& start,
                    bool triang) {
  // start position of diagonal blocks for blocked dense formats
  start.assign(n_blocks, 0);
  for (int i = 1; i < n_blocks; ++i) {
    start[i] = start[i - 1] + nb * (n - (i - 1) * nb);
    if (triang) start[i] -= nb * (nb - 1) / 2;
  }

  int jb = std::min(nb, k - (n_blocks - 1) * nb);
  double result = (double)start.back() + (double)(n - (n_blocks - 1) * nb) * jb;
  if (triang) result -= (double)jb * (jb - 1) / 2;
  return result;
}

void permuteWithSwaps(double* x, const int* swaps, int n, bool reverse) {
  // Apply swaps to vector x of length n

  if (!reverse) {
    // apply the swaps in forward order
    for (int i = 0; i < n; ++i) {
      if (swaps[i] != i) std::swap(x[i], x[swaps[i]]);
    }
  } else {
    // apply the swaps in backward order
    for (int i = n - 1; i >= 0; --i) {
      if (swaps[i] != i) std::swap(x[i], x[swaps[i]]);
    }
  }
}

void swapCols(char uplo, int n, double* A, int lda, int i, int j, int* swaps,
              int* sign) {
  // Exchange rows/cols i and j of symmetric matrix A

  // make sure that i < j
  if (i == j) return;
  if (i > j) std::swap(i, j);

  // swap diagonal elements
  std::swap(A[i + i * lda], A[j + j * lda]);

  // swap rest of rows/cols
  if (uplo == 'L') {
    callAndTime_dswap(i, &A[i], lda, &A[j], lda);
    callAndTime_dswap(n - j - 1, &A[j + 1 + i * lda], 1, &A[j + 1 + j * lda],
                      1);
    callAndTime_dswap(j - i - 1, &A[i + 1 + i * lda], 1, &A[j + (i + 1) * lda],
                      lda);
  } else {
    callAndTime_dswap(i, &A[i * lda], 1, &A[j * lda], 1);
    callAndTime_dswap(n - j - 1, &A[i + (j + 1) * lda], lda,
                      &A[j + (j + 1) * lda], lda);
    callAndTime_dswap(j - i - 1, &A[i + (i + 1) * lda], lda,
                      &A[i + 1 + j * lda], 1);
  }

  // swap pivot sign
  std::swap(sign[i], sign[j]);

  // keep track of order of swaps
  swaps[i] = j;

  DataCollector::get()->countSwap();
}

void applySwaps(const int* swaps, int nrow, int ncol, double* R) {
  // apply the column swaps to block R
  for (int i = 0; i < ncol; ++i) {
    if (swaps[i] != i) {
      // swap col i and col swaps[i]
      callAndTime_dswap(nrow, &R[i], ncol, &R[swaps[i]], ncol);
    }
  }
}

Clock::Clock() { start(); }
void Clock::start() { t0 = std::chrono::high_resolution_clock::now(); }
double Clock::stop() const {
  auto t1 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> d = t1 - t0;
  return d.count();
}
