#ifndef SENS_DOT_DRI
#define SENS_DOT_DRI

#include "../../thirdparty/asl/solvers/asl.h"
/*
 dgemv(character TRANS, integer M, integer N, double precision ALPHA, double precision, dimension(lda,*) 	A, integer LDA, double precision, dimension(*) 	X)
 integer INCX, double precision BETA, double precision, dimension(*) Y,integer INCY
*/
extern void dgemv_(char *TRANS, fint *M, fint *N, real *ALPHA, real *A, fint *LDA, real *X, fint *INCX, real *BETA, real *Y, fint *INCY);

void sens_update_driver_dot(fint n_k, fint n_u, real *s_hat_t, real *npdp, real *u_star);

#endif /* SENS_UPDT_DRI */
