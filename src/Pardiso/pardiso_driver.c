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
	fint n_rhs, real *b, real *x, fint nvar, fint ncon, int no_inertia, int nza, double logmu0){
	/* do something */
	
	int nrhs=n_rhs;
	int mtype = -2;
	void	*pt[64];

	int			iparm[64];
	double	dparm[64];
	/* real *btemp, *xtemp; */

	int maxfct, mnum, phase, error, msglvl, solver;
	int num_proc;
	int currPivotPert = 8;
	char JacPerturbed = 0;

	char *var;

	double	ddum;
	int	idum;
	int pi, ni, zi;

	double dlast, d, dmin, dmax, d0, Delta_d; 
	double km, kp, kbp;
	double kc, dc, dcb;
	int i, j, k, try_fact;
      int reduce_pivtol=0;

	double normb, normr;
	double *y;

	double nrm_r, nrm_x, nrm_b, nrm_y, ratiorr, ratiorc ;
	int incx=1;

	double const residual_ratio_max = 1e-10;
	FILE *somefile;

	dmin = 10e-20;
	dmax = 10e+40;
	d0 = 1e-08;
	kbp = 100;
	kp = 8;
	km = 1/3;
	dlast = 0.0;
	d = 0.0;
	dc = 0.0;
	try_fact = 0;

	dcb = 1e-08;
	kc = 0.25;



	error = 0;
	solver = 0;
	/* license check */
	printf("logmu0 %.g\n", logmu0);
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

	/*pardiso_printstats(&mtype, &n, a, ia, ja, &nrhs, b, &error);
	if(error != 0){
		printf("Error...\t[PRINT_STATS] %d\n", error);
		exit(1);
	}*/
	phase = 11; /* Symbolic factorization */

	iparm[10-1] = currPivotPert;
  iparm[11-1] = 2;

  iparm[13-1] = 2;

  iparm[21-1] = 3;
  iparm[24-1] = 1;
  iparm[25-1] = 1;
  iparm[29-1] = 0;
  iparm[30-1] = 1;

  dparm[ 1-1] = 300; 
  dparm[ 2-1] = 1e-6;
  dparm[ 3-1] = 5000;
  dparm[ 4-1] = 10000;
  dparm[ 5-1] = 0.5;
  dparm[ 6-1] = 0.1;
  dparm[ 7-1] = 1000;
  dparm[ 8-1] = 500;
  dparm[ 9-1] = 25;

  /*somefile = fopen("deb0.txt", "w");
  for(i=0; i<nza; i++){fprintf(somefile, "%d\t%.g\n", ja[i], a[i] );}
  fclose(somefile);
	somefile = fopen("deb01.txt", "w");
	for(i=0; i<n; i++){fprintf(somefile, "%d\n", ia[i]);}
  fclose(somefile);*/
  /* */
	pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs,
		iparm, &msglvl, &ddum, &ddum, &error, dparm);

	/*somefile = fopen("deb1.txt", "w");
  for(i=0; i<nza; i++){fprintf(somefile, "%d\t%.g\n", ja[i], a[i] );}
  fclose(somefile);
	somefile = fopen("deb11.txt", "w");
	for(i=0; i<n; i++){fprintf(somefile, "%d\n", ia[i]);}
  fclose(somefile);*/

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
  Delta_d = 0.0;
  /* check inertia
  do at least ten tries */
  for(i=0; i<10; i++){
  	phase = 11;
  	pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs,
		iparm, &msglvl, &ddum, &ddum, &error, dparm);
		phase = 22; /* Numerical factorization */
	  pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum,	&nrhs,
	  	iparm, &msglvl, &ddum, &ddum, &error, dparm);
	  /* inertias */
	  pi = iparm[21];
	  ni = iparm[22];
	  zi = n - iparm[21] - iparm[22];

	  if(error != 0){
	  	printf("Error...\t[NUM_FACT] %d\n", error);
	  	exit(2);
	  	printf("I[K_AUG]...\t[PARDISO_DRIVER]"
	  	"Factorization unsuccessful.\n");
	  }

	  

	  printf("I[K_AUG]...\t[PARDISO_DRIVER]"
	  	"Inertia (p, n, 0): (%d, %d, %d).\n", pi, ni, zi);
	  if(no_inertia){break;}

	  if(ni > ncon){
	  	fprintf(stderr, "W[K_AUG]...\t[PARDISO_DRIVER]"
	  	"Inertia check failure(neig > m).\n");
		  if(dlast == 0.0 && try_fact == 0){d = d0;}
	  	else if (dlast == 0.0 && try_fact > 0){d = kbp*d;}
	  	else if (dlast > 0.0 && try_fact > 0){d = kp*d;}
			else{d = (dmin > km*dlast) ? dmin:km*dlast;}

			if(d > dmax){exit(-1);}
			Delta_d = d - Delta_d;
			for(j=0; j<nvar; j++){
				k = ia[j]-1;
				a[k] += Delta_d; /* Add the difference */
			}
	  }
	  else if(ni < ncon){
	  	fprintf(stderr, "W[K_AUG]...\t[PARDISO_DRIVER]"
	  	"Inertia check failure(neig < m).\n");	  	
	  	if(JacPerturbed == 0){
	  		JacPerturbed = 1;
	  		fprintf(stderr, "W[K_AUG]...\t[PARDISO_DRIVER]"
	  			"Attempting to make dc > 0.\n");
	  		dc = dcb*pow((pow(10, logmu0)),kc);
	  		for(j=nvar; j<nvar+ncon; j++){
					k = ia[j]-1;
					a[k] += -dc;
				}
	  	}
	  	else{
	  		currPivotPert++;
	  		iparm[10-1] = currPivotPert;
	  	}
	  	if(currPivotPert > 12){
	  		fprintf(stderr,"E[K_AUG]...\t[PARDISO_DRIVER]"
		  		"Failure, pivot perturbation is at it maximum.\n");
	  		exit(-1);
	  	}

	  	if(JacPerturbed == 1){
	  		fprintf(stderr, "W[K_AUG]...\t[PARDISO_DRIVER]"
	  			"dc is already > 0\nThe Jacobian might be singular.\n");}
	  }
	  else if(zi!=0){
	  	fprintf(stderr,"E[K_AUG]...\t[PARDISO_DRIVER]"
		  "Failure, there is a zero eigenvalue. The jacobian is possibly singular.\n");
		  exit(-1);
		}
	  else{
	  	dlast = d;
	  	fprintf(stderr, "I[K_AUG]...\t[PARDISO_DRIVER]"
	  	"Inertia check successful.\n");
	  	/*if(try_fact > 0){break;}*/
	  	break;
  	}
		try_fact++;
  }
  phase = 33;

  printf("I[K_AUG]...\t[PARDISO_DRIVER]"
		  "Reg tries %d, reg value %.g.\n", try_fact, d);

  iparm[7] = 5;       /* Max numbers of iterative refinement steps. */ 
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
 	
 	nrm_y = dnrm2_(&n, y, &incx);
 	nrm_x = dnrm2_(&n, x, &incx);
	nrm_b = dnrm2_(&n, b, &incx);
      nrm_r = fabs(nrm_y - nrm_b);

 	/*printf("normr %e, normx %e\n", nrm_r, nrm_x);*/
 	printf("I[K_AUG]...\t[PARDISO_DRIVER]"
 		": Eucl Norm of the residuals (reported); %e \n", normr);

      
  printf("I[K_AUG]...\t[PARDISO_DRIVER]"
  	": Eucl Norm of the rhs (reported); %e \n", normb);


 	printf("I[K_AUG]...\t[PARDISO_DRIVER]"
 		": Eucl Norm of the residuals; %e \n", nrm_r);

 	printf("I[K_AUG]...\t[PARDISO_DRIVER]"
 		": Eucl Norm of the solution; %e \n", nrm_x);
      
  printf("I[K_AUG]...\t[PARDISO_DRIVER]"
  	": Eucl Norm of the rhs; %e \n", nrm_b);

      
  ratiorr = (normr/(normr + normb)) ;
  ratiorc = (nrm_r / (nrm_r + nrm_b));

	printf("I[K_AUG]...\t[PARDISO_DRIVER]"
				": Ratio of norm of scaled residuals (reported); %e \n", ratiorr);
	printf("I[K_AUG]...\t[PARDISO_DRIVER]"
				": Ratio of norm of scaled residuals (computed); %e \n", ratiorc);



	if(ratiorr > residual_ratio_max){
		printf("I[K_AUG]...\t[PARDISO_DRIVER]"
			": The norm of residuals is larger than max ratio(computed)\n");
	reduce_pivtol = 1;
      }
     	if(ratiorc > residual_ratio_max){
			printf("I[K_AUG]...\t[PARDISO_DRIVER]"
				": The norm of residuals is larger than max ratio(reported)\n");
      reduce_pivtol = 1;
	}

      if(reduce_pivtol == 1){
      		printf("I[K_AUG]...\t[PARDISO_DRIVER]"
				": Attempting to reduce residuals\n");
                  phase = 22; /* Numerical factorization */
                  pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum,	&nrhs,
            	  	iparm, &msglvl, &ddum, &ddum, &error, dparm);
                  if(error != 0){
            	  	printf("Error...\t[NUM_FACT] %d\n", error);
	  	            exit(2);
                  }

                  printf("I[K_AUG]...\t[PARDISO_DRIVER]"
            	  	"Factorization successful(new).\n");



      }


 	free(y);
	phase = -1; /* Internal memory release */
    
  pardiso (pt, &maxfct, &mnum, &mtype, &phase, &n, &ddum, ia, ja, 
  	&idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error,  dparm);

  return 0;
}
