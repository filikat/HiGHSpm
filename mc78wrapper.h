#ifndef MC78_WRAPPER_H
#define MC78_WRAPPER_H

#include "hsl_mc78i.h"

#ifdef __cplusplus
extern "C" {
#endif

void mc78_init(struct mc78_control_i* control);

int mc78_analyse(int n, const int* ptr, const int* row, int* perm, int* nnodes,
                 int** sptr, int** sparent, int64_t** rptr, int** rlist,
                 const struct mc78_control_i* control, int* stat,
                 int64_t* nfact, int64_t* nflops, int* piv_size);

#ifdef __cplusplus
}
#endif

#endif