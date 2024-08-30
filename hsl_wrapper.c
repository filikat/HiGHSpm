#include "hsl_wrapper.h"

#include <stdlib.h>

// MA57
void wrapper_ma57_default_control(struct ma57_control* control) {
  ma57_default_control(control);
}
void wrapper_ma57_init_factors(void** factors) { ma57_init_factors(factors); }
void wrapper_ma57_analyse(int n, int ne, int* row, int* col, void** factors,
                          struct ma57_control* control,
                          struct ma57_ainfo* ainfo, int* perm) {
  ma57_analyse(n, ne, row, col, factors, control, ainfo, perm);
}
void wrapper_ma57_factorize(int n, int ne, int* row, int* col, double* val,
                            void** factors, struct ma57_control* control,
                            struct ma57_finfo* finfo) {
  ma57_factorize(n, ne, row, col, val, factors, control, finfo);
}
void wrapper_ma57_solve(int n, int ne, int* row, int* col, double* val,
                        void** factors, double* x, struct ma57_control* control,
                        struct ma57_sinfo* sinfo) {
  ma57_solve(n, ne, row, col, val, factors, 1, x, control, sinfo, NULL, 0, 0);
}
void wrapper_ma57_finalise(void** factors, struct ma57_control* control,
                           int* info) {
  ma57_finalize(factors, control, info);
}

// MA86
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

// MA87
void wrapper_ma87_default_control(struct ma87_control* control) {
  ma87_default_control(control);
}

void wrapper_ma87_analyse(int n, int* ptr, int* row, int* order, void** keep,
                          struct ma87_control* control,
                          struct ma87_info* info) {
  ma87_analyse(n, ptr, row, order, keep, control, info);
}

void wrapper_ma87_factor(int n, int* ptr, int* row, double* val, int* order,
                         void** keep, struct ma87_control* control,
                         struct ma87_info* info) {
  ma87_factor(n, ptr, row, val, order, keep, control, info);
}

void wrapper_ma87_solve(int job, int nrhs, int n, double* x, int* order,
                        void** keep, struct ma87_control* control,
                        struct ma87_info* info) {
  ma87_solve(job, nrhs, n, x, order, keep, control, info);
}

void wrapper_ma87_finalise(void** keep, struct ma87_control* control) {
  ma87_finalise(keep, control);
}

// MA97
void wrapper_ma97_default_control(struct ma97_control* control) {
  ma97_default_control(control);
}

void wrapper_ma97_analyse(int n, int* ptr, int* row, int* order, void** akeep,
                          struct ma97_control* control,
                          struct ma97_info* info) {
  ma97_analyse(1, n, ptr, row, NULL, akeep, control, info, order);
}

void wrapper_ma97_factor(int type, int n, int* ptr, int* row, double* val,
                         void** akeep, void** fkeep,
                         struct ma97_control* control, struct ma97_info* info) {
  ma97_factor(type, ptr, row, val, akeep, fkeep, control, info, NULL);
}

void wrapper_ma97_solve(int job, int nrhs, int n, double* x, void** akeep,
                        void** fkeep, struct ma97_control* control,
                        struct ma97_info* info) {
  ma97_solve(job, nrhs, x, n, akeep, fkeep, control, info);
}

void wrapper_ma97_finalise(void** akeep, void** fkeep) {
  ma97_finalise(akeep, fkeep);
}

// MC68
void wrapper_mc68_default_control(struct mc68_control* control) {
  mc68_default_control(control);
}

void wrapper_mc68_order(int ord, int n, int* ptr, int* row, int* perm,
                        struct mc68_control* control, struct mc68_info* info) {
  mc68_order(ord, n, ptr, row, perm, control, info);
}