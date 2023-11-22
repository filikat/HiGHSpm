#ifndef SELECT_METHOD_H
#define SELECT_METHOD_H

#include "Metis_caller.h"
#include "util/HighsSparseMatrix.h"

int selectMethod(const HighsSparseMatrix& A, Metis_caller& Metis_data,
                  int& option_nla, int& option_metis);

int choose_NE_AS();

#endif