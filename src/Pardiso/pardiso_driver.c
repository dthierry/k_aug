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

int pardiso_driver(fint *ia, fint *ja, real *a, fint n, 
	fint n_rhs, real *b, real *x, fint nvar, fint ncon){
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
	int pi, ni, zi;

	double dlast, d, dmin, dmax, d0; 
	double km, kp, kbp;
	int i, j, k, try_fact;

	double normb, normr;
	double *y;

	double nrm_r, nrm_x;
	int incx=1;

	double const residual_ratio_max = 1e-10;

	dmin = 10e-20;
	dmax = 10e+40;
	d0 = 1e-08;
	kbp = 100;
	kp = 8;
	km = 1/3;
	dlast = 0.0;
	d = 0.0;
	try_fact = 0;


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
			printf("E[K_AUG]...\t[PARDISO_DRIVER]"
			"RHS CHECKING %d\n", error);
			exit(1);
		}

		else{
			printf("I[K_AUG]...\t[PARDISO_DRIVER]"
		"RHS CHECKING successful\n");}
	
	/* check RHS */

	pardiso_printstats(&mtype, &n, a, ia, ja, &nrhs, b, &error);
	if(error != 0){
		printf("Error...\t[PRINT_STATS] %d\n", error);
		exit(1);
	}
	phase = 11; /* Symbolic factorization */

	iparm[9] = 12;
  iparm[10] = 2;

  iparm[12] = 2;

  iparm[20] = 3;
  iparm[23] = 1;
  iparm[24] = 1;
  iparm[28] = 0;
  iparm[29] = 1;

  dparm[ 0] = 300; 
  dparm[ 1] = 1e-6;
  dparm[ 2] = 5000;
  dparm[ 3] = 10000;
  dparm[ 4] = 0.5;
  dparm[ 5] = 0.1;
  dparm[ 6] = 1000;
  dparm[ 7] = 500;
  dparm[ 8] = 25;


	pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs,
		iparm, &msglvl, &ddum, &ddum, &error, dparm);

	if (error != 0){
		printf("Error...\t[SYM_FACT] %d\n", error);
		exit(1);
	}

	printf("I[K_AUG]...\t[PARDISO_DRIVER]"
		"Reordering completed ... \n");
  printf("I[K_AUG]...\t[PARDISO_DRIVER]"
  	"Number of nonzeros in factors  = %d\n", iparm[17]);
  printf("I[K_AUG]...\t[PARDISO_DRIVER]"
  	"Number of factorization MFLOPS = %d\n", iparm[18]);

  phase = 22; /* Numerical factorization */

  /* check inertia
  do at least ten tries */
  for(i=0; i<10; i++){

	  pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum,	&nrhs,
	  	iparm, &msglvl, &ddum, &ddum, &error, dparm);
	  /* inertias */
	  pi = iparm[21];
	  ni = iparm[22];
	  zi = n - iparm[21] - iparm[22];

	  if(error != 0){
	  	printf("Error...\t[NUM_FACT] %d\n", error);
	  	exit(2);
	  }

	  printf("I[K_AUG]...\t[PARDISO_DRIVER]"
	  	"Factorization successful.\n");

	  printf("I[K_AUG]...\t[PARDISO_DRIVER]"
	  	"Inertia (p, n, 0): (%d, %d, %d).\n", pi, ni, zi);

	  if(pi != nvar || ni != ncon || zi != 0){
	  	fprintf(stderr, "W[K_AUG]...\t[PARDISO_DRIVER]"
	  	"Inertia check failure.\n");
	  }
	  else if(zi!=0){
	  	fprintf(stderr,"E[K_AUG]...\t[PARDISO_DRIVER]"
		  "Failure, there is a zero eigenvalue. The jacobian is possibly singular.\n");
		  exit(-1);
		}
	  else{
	  	dlast = d;
	  	if(try_fact > 0){break;}
  	}

  	if(dlast == 0.0 && try_fact == 0){d = d0;}
  	else if (dlast == 0.0 && try_fact > 0){d = kbp*d;}
  	else if (dlast > 0.0 && try_fact > 0){d = kp*d;}
		else{d = (dmin > km*dlast) ? dmin:km*dlast;}

		if(d > dmax){exit(-1);}

		for(j=0; j<nvar; j++){
			k = ia[j]-1;
			a[k] += d;
		}
		try_fact++;
  }
  phase = 33;

  printf("I[K_AUG]...\t[PARDISO_DRIVER]"
		  "Reg tries %d, reg value %f.\n", try_fact, d);

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
 	printf("I[K_AUG]...\t[PARDISO_DRIVER]"
 		"Solve completed x\n\n");

 	y = (double *)malloc(sizeof(double) * n * nrhs);

 	pardiso_residual (&mtype, &n, a, ia, ja, b, x, y, &normb, &normr);

 	/*printf("normb %e, normr %e\n", normb, normr);*/
 	
 	nrm_r = dnrm2_(&n, y, &incx);
 	nrm_x = dnrm2_(&n, x, &incx);
	
 	/*printf("normr %e, normx %e\n", nrm_r, nrm_x);*/


 	printf("I[K_AUG]...\t[PARDISO_DRIVER]"
			": Eucl Norm of the residuals; %e \n", nrm_r);


	printf("I[K_AUG]...\t[PARDISO_DRIVER]"
				": Ratio of norm of scaled residuals; %g \n", (nrm_r/(nrm_x+nrm_r)));
	if((nrm_r/(nrm_x+nrm_r)) < residual_ratio_max){
			printf("I[K_AUG]...\t[PARDISO_DRIVER]"
				": The norm of residuals is less than max ratio\n");
	}

 	free(y);
	phase = -1; /* Internal memory release */
    
  pardiso (pt, &maxfct, &mnum, &mtype, &phase, &n, &ddum, ia, ja, 
  	&idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error,  dparm);

  return 0;
}
