
#include "sens_update_driver_dot.h"
#include <stdlib.h>
#include <stdio.h>

/*
 dgemv(character TRANS, integer M, integer N, double precision ALPHA, double precision, dimension(lda,*) 	A, integer LDA, double precision, dimension(*) 	X)
 integer INCX, double precision BETA, double precision, dimension(*) Y,integer INCY
*/
/* int n_rhs, fint rhs_len, real *rhsbksolv, real *dp, fint nvar, fint ncon, SufDesc **rhs_ptr */
void sens_update_driver_dot(fint n_k, fint n_u, real *s_hat_t, real *npdp, real *u_star){
	char t = 'T';

	double ALPHA = -1.0;
	/*int LDA = 10;*/
	int INCX = 1;
	double BETA = 1.0;
	int INCY = 1;
  /*     y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y */
	/* Perform u = u0 - {S}^T *  Np * (p-p0) */
	/* Perform u = u0 - {S}^T *  npdp */
	dgemv_(&t, &n_k, &n_u, &ALPHA, s_hat_t, &n_k, npdp, &INCX, &BETA, u_star, &INCY);
	/*dgemv_(&t, &n_u, &n_k, &ALPHA, s_hat_t, &n_u, npdp, &INCX, &BETA, u_star, &INCY);	*/

}