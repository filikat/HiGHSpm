#ifndef VERTEX_COVER_H
#define VERTEX_COVER_H

#include <vector>
#include <set>
#include <random>
#include <algorithm>
#include <climits>

void vertexCoverMM(int nvertex, int nedges, int nparts,
                   const std::vector<int> part, const std::vector<int> adj_ptr,
                   const std::vector<int> adj_lst,
                   std::vector<int>& permutation, std::vector<int>& blockSize);

void vertexCoverG(int nvertex, int nedges, int nparts,
                  const std::vector<int> part, const std::vector<int> adj_ptr,
                  const std::vector<int> adj_lst, std::vector<int>& permutation,
                  std::vector<int>& blockSize);

#endif