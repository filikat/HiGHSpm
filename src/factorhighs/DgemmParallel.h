#ifndef DGEMM_PARALLEL_H
#define DGEMM_PARALLEL_H

// parallelize dgemm for use within factorization
// Performs Q <- Q - R P^T in hybrid format.
// Parallelized over the rows of R and Q.
class dgemmParallelizer {
  const double* P_;
  const double* R_;
  double* Q_;
  const int col_;
  const int jb_;

 public:
  dgemmParallelizer(const double* P, const double* R, double* Q, int col,
                    int jb);

  void run(int start, int end, double beta) const;
};

void dgemmParallel(const double* P, const double* R, double* Q, int col, int jb,
                   int row, int nb, double beta = 1.0);

#endif