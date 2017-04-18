
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "asl.h"

// Prototypes
void pardisoinit (void   *, int    *,   int *, int *, double *, int *);
void pardiso     (void   *, int    *,   int *, int *,    int *, int *, 
                  double *, int    *,    int *, int *,   int *, int *,
                     int *, double *, double *, int *, double *);
void pardiso_chkmatrix  (int *, int *, double *, int *, int *, int *);
void pardiso_chkvec     (int *, int *, double *, int *);
void pardiso_printstats (int *, int *, double *, int *, int *, int *,
                           double *, int *);

int pardiso_driver(fint *ia, fint *ja, real *a, fint n, fint nza,
 fint nrhs, real *b, real *x);

int pardiso_driver(fint *ia, fint *ja, real *a, fint n, fint nza, 
	fint nrhs, real *b, real *x){
	// do something
	int mtype = -2;
	void	*pt[64];

	int			iparm[64];
	double	dparm[64];

	int maxfct, mnum, phase, error, msglvl, solver;
	int num_proc;

	char *var;
	int i, k;

	double	ddum;
	int			idum;

	error = 0;
	solver = 0;

	pardisoinit(pt, &mtype, &solver, iparm, dparm, &error);

	if(error != 0){
		if(error == -10)
			printf("No license file found \n");
		else if (error == -11)
			printf("License is expired \n");
		else if (error = -12)
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

	msglvl = 1;
	error = 0;

	pardiso_chkmatrix(&mtype, &n, a, ia, ja, &error);
	if(error != 0){
		printf("Error in...\t[MATRIX] %d\n", error);
		exit(1);
	}


	pardiso_chkvec(&n, &nrhs, b, &error);

	if(error != 0){
		printf("Error...\t[RHS] %d\n", error);
		exit(1);
	}

	pardiso_printstats(&mtype, &n, a, ia, ja, &nrhs, b, &error);
	if(error != 0){
		printf("Error...\t[PRINT_STATS] %d\n", error);
		exit(1);
	}

	phase = 11;

	pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs,
		iparm, &msglvl, &ddum, &ddum, &error, dparm);

	if (error != 0){
		printf("Error...\t[SYM_FACT] %d\n", error);
		exit(1);
	}

	printf("Reordering completed ... \n");
  printf("Number of nonzeros in factors  = %d\n", iparm[17]);
  printf("\nNumber of factorization MFLOPS = %d\n", iparm[18]);

  phase = 22;
  pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum,	&nrhs,
  	iparm, &msglvl, &ddum, &ddum, &error, dparm);

  if(error != 0){
  	printf("Error...\t[NUM_FACT] %d\n", error);
  	exit(2);
  }

  printf("\nFactorization successful!\n");

	phase = -1;                 /* Release internal memory. */
    
  pardiso (pt, &maxfct, &mnum, &mtype, &phase, &n, &ddum, ia, ja, 
  	&idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error,  dparm);

  return 0;
}