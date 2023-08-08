#ifndef HSL_MA86_WRAPPER_H
#define HSL_MA86_WRAPPER_H
#include "hsl_mc69d.h"
#include "hsl_ma86d.h"

#ifdef __cplusplus
extern "C" {
#endif

void wrapper_ma86_default_control(struct ma86_control_d* control);
void wrapper_ma86_analyse(int n, int* ptr, int* row, int* order, void** keep, struct ma86_control_d* control, struct ma86_info_d* info);
void wrapper_ma86_factor(int n, int* ptr, int* row, double* val, int* order, void** keep, struct ma86_control_d* control, struct ma86_info_d* info);
void wrapper_ma86_solve(int job, int nrhs, int n, double* x, int* order, void** keep, struct ma86_control_d* control, struct ma86_info_d* info);
void wrapper_ma86_finalise(void** keep, struct ma86_control_d* control);

#ifdef __cplusplus
}
#endif

#endif  // HSL_MA86_WRAPPER_H
