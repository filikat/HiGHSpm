#include "SparseMatrix.h"


    // =======================================================================
    // CONSTRUCTOR
    // =======================================================================
    SparseMatrix::SparseMatrix(int num_row,int num_col,int num_nnz){
        row_index.resize(num_nnz);
        col_ptr.resize(num_col+1);
        values.resize(num_nnz);
        m = num_row;
        n = num_col;
    }



// =======================================================================
// MATRIX-VECTOR PRODUCT 
// =======================================================================
void mat_vec( 
    const SparseMatrix& A,
    const std::vector<double>& x,
    std::vector<double>& y,
    double alpha,
    char tran = 'n'
    )
{
    if (tran == 't' || tran == 'T'){
        // loop through columns of A
        for (int i=0; i<A.cols(); ++i){
            //loop through nonzero entries of column
            for (int j=A.begin(i); j<A.end(i); ++j){
                y[i] += alpha * A.values[j] * x[A.row_index[j]];
            }
        }
    }else{
        // loop through columns of A
        for (int i=0; i<A.cols(); ++i){
            //loop through nonzero entries of column
            for (int j=A.begin(i); j<A.end(i); ++j){
                y[A.row_index[j]] += alpha * A.values[j] * x[i];
            }
        }
    }
}


// =======================================================================
// TRANSPOSE
// =======================================================================
SparseMatrix SparseMatrix::transpose() const{

    SparseMatrix T(n,m,nnz());

    // count elements in each row
    std::vector<int> work(m,0);
    for (int i=0; i<nnz(); ++i){
        ++work[row_index[i]];
    }

    // accumulate rows and create row pointers
    T.col_ptr[0] = 0;
    for (int i=0; i<m; ++i){
        T.col_ptr[i+1] = T.col_ptr[i] + work[i];

        // w contains starting point of each row and will be incremented
        work[i] = T.col_ptr[i];
    }

    //go through the columns
    for (int i=0; i<n; ++i){
        for (int j=begin(i); j<end(i); ++j){

            // work[row_index[j]] is the next available position for entries in column row_index[j] of the transpose matrix
            T.row_index[work[row_index[j]]] = i;
            T.values[work[row_index[j]]] = values[j];
            ++work[row_index[j]];
        }
    }

    return T;
}



std::ostream& operator<<(std::ostream& out,const SparseMatrix& M){

    out << "row indices: ";
    for (int i: M.row_index){
        out << i << ' ';
    }
    out << "\ncol ptr:     ";
    for (int i: M.col_ptr){
        out << i << ' ';
    }
    out << "\nvalues:      ";
    for (double i: M.values){
        out << i << ' ';
    }

    return out;
}



