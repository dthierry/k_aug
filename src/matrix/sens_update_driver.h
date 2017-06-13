#ifndef SENS_UPDT_DRI
#define SENS_UPDT_DRI

#include "asl.h"
/*
 dgemv(character TRANS, integer M, integer N, double precision ALPHA, double precision, dimension(lda,*) 	A, integer LDA, double precision, dimension(*) 	X)
 integer INCX, double precision BETA, double precision, dimension(*) Y,integer INCY
*/
extern void dgemv_(char *TRANS, int *M, int *N, double *ALPHA, double *A, int *LDA, double *X, int *INCX, double *BETA, double *Y, int *INCY);

void sens_update_driver(int n_rhs, fint rhs_len, real *dsdp, real *dp, real *sstar);

#endif /* SENS_UPDT_DRI */