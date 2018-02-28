/* @source test_mumps.c
**
** January 12th, 2018
** @author: David Thierry (dmolinat@andrew.cmu) dav0@lb2016-1

******************************************************************************

@ptest_mumps ********************************************
**
** 
**


** @@ This comes mostly from the example file at their website
*******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include "mumps_driver.h"


int main (int argc, char **argc){
	int n=5;
	int nnz = 12;
	int irn[] = {1, 2, 4, 5, 2, 1, 5, 3, 2, 3, 1, 3};
	int jcn[] = {2, 3, 3, 5, 1, 1, 2, 4, 5, 2, 3, 3};
	double a[] = {3.0, -3.0, 2.0, 1.0, 3.0, 2.0, 4.0, 2.0, 6.0, -1.0, 4.0, 1.0};
	double b[] = {20.0, 24.0, 9.0, 6.0, 13.0};
	int n_rhs = 1;
	double *x=NULL;
	int nvar=5, ncon=5;
	int no_inertia=1;
	double logmu0 = -8.99;

	mumps_driver(irn, jcn, a, &n, &n_rhs, b, x, nvar, ncon, no_inertia, &nnz, &logmu0);

	return 0;
}