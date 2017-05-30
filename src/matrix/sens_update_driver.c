
#include "sens_update_driver.h"
#include <stdlib.h>
#include <stdio.h>

/*
 dgemv(character TRANS, integer M, integer N, double precision ALPHA, double precision, dimension(lda,*) 	A, integer LDA, double precision, dimension(*) 	X)
 integer INCX, double precision BETA, double precision, dimension(*) Y,integer INCY
*/
/* int n_rhs, fint rhs_len, real *rhsbksolv, real *dp, fint nvar, fint ncon, SufDesc **rhs_ptr */
void sens_update_driver(int n_rhs, fint rhs_len, real *dsdp, real *dp, real *sstar){
	char t = 'N';

	double ALPHA = 1.0;
	/*int LDA = 10;*/
	int INCX = 1;
	double BETA = 1.0;
	int INCY = 1;


	dgemv_(&t, &rhs_len, &n_rhs, &ALPHA, dsdp, &rhs_len, dp, &INCX, &BETA, sstar, &INCY);	

}