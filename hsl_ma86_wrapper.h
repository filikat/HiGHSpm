#ifndef HSL_MA86_WRAPPER_H
#define HSL_MA86_WRAPPER_H
#include "hsl_ma86d.h"
#include "hsl_mc68i.h"
#include "hsl_mc69d.h"

#ifdef __cplusplus
extern "C" {
#endif

void wrapper_ma86_default_control(struct ma86_control* control);
void wrapper_ma86_analyse(int n, int* ptr, int* row, int* order, void** keep,
                          struct ma86_control* control, struct ma86_info* info);
void wrapper_ma86_factor(int n, int* ptr, int* row, double* val, int* order,
                         void** keep, struct ma86_control* control,
                         struct ma86_info* info);
void wrapper_ma86_solve(int job, int nrhs, int n, double* x, int* order,
                        void** keep, struct ma86_control* control,
                        struct ma86_info* info);
void wrapper_ma86_finalise(void** keep, struct ma86_control* control);

void wrapper_mc68_default_control(struct mc68_control* control);
void wrapper_mc68_order(int ord, int n, int* ptr, int* row, int* perm,
                        struct mc68_control* control, struct mc68_info* info);

#ifdef __cplusplus
}
#endif

#endif  // HSL_MA86_WRAPPER_H
