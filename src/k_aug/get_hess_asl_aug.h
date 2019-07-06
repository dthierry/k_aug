#ifndef HESS_AUG
#define HESS_AUG

#include "../../thirdparty/asl/solvers/asl.h"


void get_hess_asl_aug(ASL *asl, real *x, fint **Wcol, fint **Wrow, real **Wij,
                      int nvar, int ncon, int nobj, fint *n_nz_w, real *y, fint *nerror,  int **nz_row_w,
                      int **md_off_w, int *missing_nz);

#endif /* HESS_AUG */
