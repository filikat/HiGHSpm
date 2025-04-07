#ifndef AUXILIARY_H
#define AUXILIARY_H

#include <chrono>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

#include "util/HighsCDouble.h"

void counts2Ptr(std::vector<int>& ptr, std::vector<int>& w);
void inversePerm(const std::vector<int>& perm, std::vector<int>& iperm);
void subtreeSize(const std::vector<int>& parent, std::vector<int>& sizes);
void transpose(const std::vector<int>& ptr, const std::vector<int>& rows,
               std::vector<int>& ptrT, std::vector<int>& rowsT);
void transpose(const std::vector<int>& ptr, const std::vector<int>& rows,
               const std::vector<double>& val, std::vector<int>& ptrT,
               std::vector<int>& rowsT, std::vector<double>& valT);
void symProduct(const std::vector<int>& ptr, const std::vector<int>& rows,
                const std::vector<double>& vals, const std::vector<double>& x,
                std::vector<double>& y, double alpha = 1.0);
void symProductQuad(const std::vector<int>& ptr, const std::vector<int>& rows,
                    const std::vector<double>& vals,
                    const std::vector<double>& x, std::vector<HighsCDouble>& y,
                    double alpha);
void childrenLinkedList(const std::vector<int>& parent, std::vector<int>& head,
                        std::vector<int>& next);
void reverseLinkedList(std::vector<int>& head, std::vector<int>& next);
void dfsPostorder(int node, int& start, std::vector<int>& head,
                  const std::vector<int>& next, std::vector<int>& order);
void processEdge(int j, int i, const std::vector<int>& first,
                 std::vector<int>& maxfirst, std::vector<int>& delta,
                 std::vector<int>& prevleaf, std::vector<int>& ancestor);
double getDiagStart(int n, int k, int nb, int n_blocks, std::vector<int>& start,
                    bool triang = false);
void permuteWithSwaps(double* x, const int* swaps, int n, bool reverse = false);
void swapCols(char uplo, int n, double* A, int lda, int i, int j, int* swaps,
              int* sign);
void applySwaps(const int* swaps, int nrow, int ncol, double* R);

template <typename T>
void permuteVector(std::vector<T>& v, const std::vector<int>& perm) {
  // Permute vector v according to permutation perm.
  std::vector<T> new_v(v.size());
  for (int i = 0; i < v.size(); ++i) {
    new_v[i] = v[perm[i]];
  }
  v = std::move(new_v);
}

template <typename T>
void permuteVectorInverse(std::vector<T>& v, const std::vector<int>& iperm) {
  // Permute vector v according to inverse permutation iperm.
  std::vector<T> new_v(v.size());
  for (int i = 0; i < v.size(); ++i) {
    new_v[iperm[i]] = v[i];
  }
  v = std::move(new_v);
}

template <typename T>
void print(const std::vector<T>& v, const std::string s) {
  std::ofstream out_file;
  char name[80];
  snprintf(name, 80, "../FactorHiGHS/matlab/%s.txt", s.c_str());
  out_file.open(name);
  for (T i : v) {
    out_file << std::setprecision(16) << i << '\n';
  }
  out_file.close();
}

class Clock {
  std::chrono::high_resolution_clock::time_point t0;

 public:
  Clock();
  void start();
  double stop() const;
};

#endif
