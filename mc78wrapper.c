#include "mc78wrapper.h"

void mc78_init(struct mc78_control_i* control) {
  mc78_default_control(control);
}

int mc78_analyse(int n, const int* ptr, const int* row, int* perm, int* nnodes,
                 int** sptr, int** sparent, int64_t** rptr, int** rlist,
                 const struct mc78_control_i* control, int* stat,
                 int64_t* nfact, int64_t* nflops, int* piv_size) {
  return mc78_analyse_asm(n, ptr, row, perm, nnodes, sptr, sparent, rptr, rlist,
                          control, stat, nfact, nflops, piv_size);
}