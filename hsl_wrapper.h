#ifndef HSL_WRAPPER_H
#define HSL_WRAPPER_H
#include "hsl_ma57d.h"
#include "hsl_ma86d.h"
#include "hsl_ma87d.h"
#include "hsl_ma97d.h"
#include "hsl_mc68i.h"
#include "hsl_mc69d.h"

#ifdef __cplusplus
extern "C" {
#endif

// MA57
void wrapper_ma57_default_control(struct ma57_control *control);
void wrapper_ma57_init_factors(void **factors);
void wrapper_ma57_analyse(int n, int ne, int *row, int *col, void **factors,
                          struct ma57_control *control,
                          struct ma57_ainfo *ainfo, int *perm);

void wrapper_ma57_factorize(int n, int ne, int *row, int *col, double *val,
                            void **factors, struct ma57_control *control,
                            struct ma57_finfo *finfo);
void wrapper_ma57_solve(int n, int ne, int *row, int *col, double *val,
                        void **factors, double *x, struct ma57_control *control,
                        struct ma57_sinfo *sinfo);
void wrapper_ma57_finalise(void **factors, struct ma57_control *control,
                           int *info);

// MA86
void wrapper_ma86_default_control(struct ma86_control *control);
void wrapper_ma86_analyse(int n, int *ptr, int *row, int *order, void **keep,
                          struct ma86_control *control, struct ma86_info *info);
void wrapper_ma86_factor(int n, int *ptr, int *row, double *val, int *order,
                         void **keep, struct ma86_control *control,
                         struct ma86_info *info);
void wrapper_ma86_solve(int job, int nrhs, int n, double *x, int *order,
                        void **keep, struct ma86_control *control,
                        struct ma86_info *info);
void wrapper_ma86_finalise(void **keep, struct ma86_control *control);

// MA87
void wrapper_ma87_default_control(struct ma87_control *control);
void wrapper_ma87_analyse(int n, int *ptr, int *row, int *order, void **keep,
                          struct ma87_control *control, struct ma87_info *info);
void wrapper_ma87_factor(int n, int *ptr, int *row, double *val, int *order,
                         void **keep, struct ma87_control *control,
                         struct ma87_info *info);
void wrapper_ma87_solve(int job, int nrhs, int n, double *x, int *order,
                        void **keep, struct ma87_control *control,
                        struct ma87_info *info);
void wrapper_ma87_finalise(void **keep, struct ma87_control *control);

// MA97
void wrapper_ma97_default_control(struct ma97_control *control);
void wrapper_ma97_analyse(int n, int *ptr, int *row, int *order, void **akeep,
                          struct ma97_control *control, struct ma97_info *info);
void wrapper_ma97_factor(int type, int n, int *ptr, int *row, double *val,
                         void **akeep, void **fkeep,
                         struct ma97_control *control, struct ma97_info *info);
void wrapper_ma97_solve(int job, int nrhs, int n, double *x, void **akeep,
                        void **fkeep, struct ma97_control *control,
                        struct ma97_info *info);
void wrapper_ma97_finalise(void **akeep, void **fkeep);

// MA68
void wrapper_mc68_default_control(struct mc68_control *control);
void wrapper_mc68_order(int ord, int n, int *ptr, int *row, int *perm,
                        struct mc68_control *control, struct mc68_info *info);

#ifdef __cplusplus
}
#endif

#endif
