/* @source pardiso_driver.c
**
** April 18th, 2017
** @author: David Thierry (dmolinat@andrew.cmu) dav0@lb2016-1

******************************************************************************

@pardiso_driver ********************************************
**
** calls pardiso, solves linear system
**
** @param [r] ia, row starts
** @param [r] ja, sparse column
** @param [r] a, sparse matrix element
** @param [r] n, number of rows
** @param [r] nza, number of nz of a
** @param [r] nrhs, number of right hand sides
** @param [r] b, right hand side vector(matrix)
** @param [r] x, solution vector

** @@ This comes mostly from the example file at their website
*******************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "pardiso_driver.h"

int pardiso_driver(fint *ia, fint *ja, real *a, fint n, fint nza, 
	fint n_rhs, real *b, real *x){
	/* do something */
	
	int nrhs=n_rhs;
	int mtype = -2;
	void	*pt[64];

	int			iparm[64];
	double	dparm[64];
	/* real *btemp, *xtemp; */

	int maxfct, mnum, phase, error, msglvl, solver;
	int num_proc;

	char *var;

	double	ddum;
	int			idum;

	error = 0;
	solver = 0;
	/* license check */
	pardisoinit(pt, &mtype, &solver, iparm, dparm, &error);

	if(error != 0){
		if(error == -10)
			printf("No license file found \n");
		else if (error == -11)
			printf("License is expired \n");
		else if (error == -12)
			printf("Wrong username or hostname \n");
		return 1;
	}
	else
		printf("PARDISO.. \t License check was successful\n");

	num_proc = 1;
	var = getenv("OMP_NUM_THREADS");

	if(var != NULL)
		sscanf( var, "%d", &num_proc);
	else{
		printf("Set environment OMP_NUM_THREADS to 1\n");
	}

	iparm[2] = num_proc;

	maxfct = 1;
	mnum = 1;

	msglvl = 0;
	error = 0;

	/* check Matrix */
	pardiso_chkmatrix(&mtype, &n, a, ia, ja, &error);
	if(error != 0){
		printf("Error in...\t[MATRIX] %d\n", error);
		exit(1);
	}


		pardiso_chkvec(&n, &nrhs, b, &error);
		if(error != 0){
			printf("E[KMATRIX]...\t[PARDISO_DRIVER]"
			"RHS CHECKING %d\n", error);
			exit(1);
		}

		else{
			printf("I[KMATRIX]...\t[PARDISO_DRIVER]"
		"RHS CHECKING successful\n");}
	
	/* check RHS */

	pardiso_printstats(&mtype, &n, a, ia, ja, &nrhs, b, &error);
	if(error != 0){
		printf("Error...\t[PRINT_STATS] %d\n", error);
		exit(1);
	}
	phase = 11; /* Symbolic factorization */

	pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs,
		iparm, &msglvl, &ddum, &ddum, &error, dparm);

	if (error != 0){
		printf("Error...\t[SYM_FACT] %d\n", error);
		exit(1);
	}

	printf("Reordering completed ... \n");
  printf("Number of nonzeros in factors  = %d\n", iparm[17]);
  printf("\nNumber of factorization MFLOPS = %d\n", iparm[18]);

  phase = 22; /* Numerical factorization */
  pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum,	&nrhs,
  	iparm, &msglvl, &ddum, &ddum, &error, dparm);

  if(error != 0){
  	printf("Error...\t[NUM_FACT] %d\n", error);
  	exit(2);
  }

  printf("I[KMATRIX]...\t[PARDISO_DRIVER]"
  	"Factorization successful.\n");

  printf("I[KMATRIX]...\t[PARDISO_DRIVER]"
  	"Inertia (p, n, 0): (%d, %d, %d).\n", iparm[21], iparm[22], n - iparm[21] - iparm[22]);

  phase = 33;


  iparm[7] = 1;       /* Max numbers of iterative refinement steps. */ 
  /*btemp = b; */
  /*xtemp = *x; */
  /*exit(1); */
	pardiso (pt, &maxfct, &mnum, &mtype, &phase,
           &n, a, ia, ja, &idum, &nrhs,
           iparm, &msglvl, b, x, &error,  dparm);
 
  if (error != 0) {
    printf("\nERROR during solution: %d", error);
    exit(3);
  }
 	printf("I[KMATRIX]...\t[PARDISO_DRIVER]"
 		"Solve completed x\n\n");
	

	phase = -1; /* Internal memory release */
    
  pardiso (pt, &maxfct, &mnum, &mtype, &phase, &n, &ddum, ia, ja, 
  	&idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error,  dparm);

  return 0;
}
