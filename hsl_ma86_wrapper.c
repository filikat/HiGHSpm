#include "hsl_ma86_wrapper.h"

#include <stdlib.h>

#include "hsl_ma86d.h"
#include "hsl_mc69d.h"

void wrapper_ma86_default_control(struct ma86_control* control) {
  ma86_default_control(control);
}

void wrapper_ma86_analyse(int n, int* ptr, int* row, int* order, void** keep,
                          struct ma86_control* control,
                          struct ma86_info* info) {
  ma86_analyse(n, ptr, row, order, keep, control, info);
}

void wrapper_ma86_factor(int n, int* ptr, int* row, double* val, int* order,
                         void** keep, struct ma86_control* control,
                         struct ma86_info* info) {
  ma86_factor(n, ptr, row, val, order, keep, control, info, NULL);
}

void wrapper_ma86_solve(int job, int nrhs, int n, double* x, int* order,
                        void** keep, struct ma86_control* control,
                        struct ma86_info* info) {
  ma86_solve(job, nrhs, n, x, order, keep, control, info, NULL);
}

void wrapper_ma86_finalise(void** keep, struct ma86_control* control) {
  ma86_finalise(keep, control);
}

void wrapper_mc68_default_control(struct mc68_control* control) {
  mc68_default_control(control);
}

void wrapper_mc68_order(int ord, int n, int* ptr, int* row, int* perm,
                        struct mc68_control* control, struct mc68_info* info) {
  mc68_order(ord, n, ptr, row, perm, control, info);
}