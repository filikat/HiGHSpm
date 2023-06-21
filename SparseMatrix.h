#ifndef SPARSEMATRIX_H
#define SPARSEMATRIX_H

#include <iostream>
#include <vector>

class SparseMatrix {

  std::vector<int> row_index{};
  std::vector<int> col_ptr{};
  std::vector<double> values{};

  int m{}; // rows
  int n{}; // columns

public:
  // =======================================================================
  // CONSTRUCTORS
  // =======================================================================

  SparseMatrix() = default;

  // Allocate space for matrix
  SparseMatrix(int num_row, int num_col, int num_nnz);

  // copy data into matrix
  SparseMatrix(const std::vector<int> &input_row_index,
               const std::vector<int> &input_col_ptr,
               const std::vector<double> &input_values, int input_m,
               int input_n)
      : row_index{input_row_index}, col_ptr{input_col_ptr},
        values{input_values}, m{input_m}, n{input_n} {}

  // =======================================================================
  // INFO
  // =======================================================================
  int rows() const { return m; }
  int cols() const { return n; }
  int nnz() const { return col_ptr.back(); }
  int begin(int j) const { return col_ptr[j]; }
  int end(int j) const { return col_ptr[j + 1]; }
  int Row(int j) const { return row_index[j]; }
  int Val(int j) const { return values[j]; }

  // =======================================================================
  // TRANSPOSE
  // =======================================================================
  SparseMatrix transpose() const;

  // =======================================================================
  // MATRIX-VECTOR PRODUCT
  // =======================================================================
  // Compute  y += alpha * A * x
  //          OR
  //          y += alpha * A^T * x , if tran is 't' or 'T'
  //
  // =======================================================================
  friend void mat_vec(              // length
      const SparseMatrix &A,        // m x n
      const std::vector<double> &x, // n (or m if 'T')
      std::vector<double> &y,       // m (or n if 'T')
      double alpha, char tran);

  friend std::ostream &operator<<(std::ostream &, const SparseMatrix &);
};

#endif