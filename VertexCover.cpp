#include "VertexCover.h"

void vertexCoverMM(int nvertex, int nedges, int nparts,
                   const std::vector<int> part, const std::vector<int> adj_ptr,
                   const std::vector<int> adj_lst,
                   std::vector<int>& permutation, std::vector<int>& blockSize) {
  // ptr and adj of the nodes in the cut
  std::vector<int> cut_ptr;
  cut_ptr.reserve(nvertex + 1);
  cut_ptr.push_back(0);
  std::vector<int> cut_adj;
  cut_adj.reserve(nedges / 10);

  // keep track of which nodes are in the cut
  std::vector<bool> cut_nodes(nvertex, false);

  // blocks of nodes with the same partition
  std::vector<std::set<int>> blocks(nparts, std::set<int>{});

  // go through the nodes
  for (int node = 0; node < nvertex; ++node) {
    int partNode = part[node];
    blocks[partNode].insert(node);

    int countCutNeigh{};

    // go through neighbours of node
    for (int j = adj_ptr[node]; j < adj_ptr[node + 1]; ++j) {
      int neigh = adj_lst[j];

      // neighbour has different partition
      if (partNode != part[neigh]) {
        ++countCutNeigh;
        cut_adj.push_back(neigh);
        cut_nodes[node] = true;
      }
    }

    cut_ptr.push_back(cut_ptr.back() + countCutNeigh);
  }

  // randon number generator
  std::random_device rd;
  std::mt19937 rng{rd()};

  std::set<int> bestNodeCover{};
  int bestNodeCoverSize{INT_MAX};

  // perform multiple maximal matching and choose the best
  for (int iter = 0; iter < 10; ++iter) {
    // random order in which to visit the vertices
    std::vector<int> index(nvertex);
    for (int i = 0; i < nvertex; ++i) {
      index[i] = i;
    }
    std::shuffle(index.begin(), index.end(), rng);

    // shuffle adjacency of each node
    for (int i = 0; i < nvertex; ++i) {
      std::shuffle(cut_adj.begin() + cut_ptr[i],
                   cut_adj.begin() + cut_ptr[i + 1], rng);
    }

    std::set<int> nodeCut{};

    std::vector<bool> matched(nvertex, false);

    // go through the vertices
    for (int i = 0; i < nvertex; ++i) {
      // select a random node
      int node = index[i];

      // skip node if it's not in the cut
      if (!cut_nodes[node]) {
        continue;
      }

      // node is not available anymore to be matched
      matched[node] = true;

      // go through the neighbours
      for (int j = cut_ptr[node]; j < cut_ptr[node + 1]; ++j) {
        int neigh = cut_adj[j];

        // found an available neighbour to match
        if (matched[neigh] == false) {
          matched[neigh] = true;
          nodeCut.insert(node);
          nodeCut.insert(neigh);
          break;
        }
      }
    }

    if (nodeCut.size() < bestNodeCoverSize) {
      bestNodeCoverSize = nodeCut.size();
      bestNodeCover = nodeCut;
    }
  }

  // remove nodes in the cut from the permutation
  for (int i : bestNodeCover) {
    blocks[part[i]].erase(i);
  }

  permutation.resize(nvertex);
  blockSize.resize(nparts + 1);

  int current{};
  for (int i = 0; i < nparts; ++i) {
    for (int j : blocks[i]) {
      permutation[current] = j;
      ++current;
    }
    blockSize[i] = blocks[i].size();
  }
  for (int i : bestNodeCover) {
    permutation[current] = i;
    ++current;
  }
  blockSize[nparts] = bestNodeCover.size();
}

void vertexCoverN(int nvertex, int nedges, int nparts,
                  const std::vector<int> part, const std::vector<int> adj_ptr,
                  const std::vector<int> adj_lst, std::vector<int>& permutation,
                  std::vector<int>& blockSize) {
  // set up stuff
  std::vector<bool> nodes_removed(nvertex, false);
  std::vector<std::set<int>> blocks(nparts, std::set<int>{});
  std::set<int> linking{};

  // go through the nodes
  for (int node = 0; node < nvertex; ++node) {
    int partNode = part[node];

    // go through the neighbours
    for (int j = adj_ptr[node]; j < adj_ptr[node + 1]; ++j) {
      int neigh = adj_lst[j];

      // skip selp loops
      if (neigh == node) continue;

      // skip if neighbour has the same partition
      if (partNode == part[neigh]) {
        continue;
      }

      // skip if neighbour is yet to be visited
      if (neigh > node) {
        continue;
      }

      // if neither node or neighbour have been removed
      if (!nodes_removed[node]) {
        if (!nodes_removed[neigh]) {
          nodes_removed[node] = true;
          linking.insert(node);
        }
      }
    }

    blocks[part[node]].insert(node);
  }

  // remove nodes in the cut from the permutation
  for (int i : linking) {
    blocks[part[i]].erase(i);
  }

  permutation.resize(nvertex);
  blockSize.resize(nparts + 1);

  int current{};
  for (int i = 0; i < nparts; ++i) {
    for (int j : blocks[i]) {
      permutation[current] = j;
      ++current;
    }
    blockSize[i] = blocks[i].size();
  }
  for (int i : linking) {
    permutation[current] = i;
    ++current;
  }
  blockSize[nparts] = linking.size();
}