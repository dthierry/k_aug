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

#include "../../ASL/solvers/asl.h"
#include "mumps_driver.h"


int main (int argc, char **argv){
	int i;
    int n = 3; /*! something interesting here */
	int nnz = 6;
	int irn[] = {1, 1, 1, 2, 2, 3};
	int jcn[] = {1, 2, 3, 2, 3, 3};
	double a[] = {5.0, 6.0, 7.0, 3.0, 2.0, 1.0};
	double b[] = {20.0, 24.0, 9.0};
	int n_rhs = 1;
	double *x=NULL;
	int nvar=3, ncon=3;
	int no_inertia=1;
	double logmu0 = -8.99;
	int retval;

	retval = mumps_driver(irn, jcn, a, n, n_rhs, b, x, nvar, ncon, no_inertia, nnz, logmu0);
    for(i=0; i<3; i++){
        printf("result %8.2f\n", *(b + i));
    }
	if(retval==0){
        printf("Success\n");
	}
	return 0;
}