/* @source dpotri_driver.c
** beta 01
** June 26th, 2017
** @author: David Thierry (dmolinat@andrew.cmu) dav0@lb2016-1

********************************************************************************

@fun_name ********************************************
**
** Description
** Description
**
** @param [r] 
** @param [r] 
** @param [r] 
** @param [r] 
** @param [r] 
** @param [r] 
** @param [r] 
** @return something
*******************************************************************************/

#include "dpotri_driver.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

double max_fun(int n, double *_ai, int *imax);

void dpotri_driver(int n, double *_a, long Kn, int *sb_p, char *_chr_timest){
	char uplo;
	double *w, *work;
	double *a;
	double *ai; /* For the modified Cholesky */
	double amax=0.0, *aimax, xi, sqbeta, delta, ainf;
	int imax=0, dummy=0;
	double w_mock[1], work_mock[1];
	int lwork, lda;
	int info=0;
	int i, j, k, try_d;
	int ret_val=0;
	FILE *somefile;
	char _file_name_[30] = {"inv_"};
	double dlast, d, dmin, dmax, d0; 
	double km, kp, kbp;

	strcat(_file_name_, _chr_timest);
	strcat(_file_name_, ".in");
	fprintf(stderr, "I[K_AUG]...\t[DPOTRI_DRIVER] Output file name %s\n", _file_name_);

 	/* transform a to up-triangular */
	
	uplo = 'U'; /* Uppertriangular */
	lda = n;
	lwork = -1;

	dmin = 10e-20;
	dmax = 10e+12;
	d0 = 1e-03;
	kbp = 100;
	kp = 100;
	km = 1/3;
	dlast = 0.0;
	d = 0.0;
	try_d = 0;

	a = (double *)malloc(sizeof(double) * n * n);
	ai = (double *)calloc(sizeof(double), n);
	aimax = (double *)calloc(sizeof(double), n);


	for(i=0; i<n; i++){
		for(j=0; j<n; j++){
			if(j > i){
				*(a + n * i + j) = 0.0;
			}
			else{
				*(a + n * i + j) = (*(_a + Kn * i + sb_p[j]) + *(_a + Kn * j + sb_p[i])) * 0.5;
			}
			if(i == j){ai[i] = *(a + n * i + j);}
		}
	}
	somefile = fopen("inv_red_hess", "w");
	for(j=0; j<n; j++){
        for(k=0; k<n; k++){
            fprintf(somefile, "%.10g\t", *(a + n * k + j));
        }
        fprintf(somefile, "\n");}
	fclose(somefile);

	
	/* Since the inertia was previously checked I expect the eigenvalues to be positive, and, if there are
	negative ones; very small. */
	
	

	for(i=0; i<10; i++){
		somefile = fopen("a0", "w");
		for(j=0; j<n; j++){for(k=0; k<n; k++){fprintf(somefile, "%f\t", *(a + n * k + j));}fprintf(somefile, "\n");}
		fclose(somefile);
		dpotrf_(&uplo, &n, a, &lda, &info); /* Attempt once */
		somefile = fopen("a1", "w");
		for(j=0; j<n; j++){for(k=0; k<n; k++){fprintf(somefile, "%f\t", *(a + n * k + j));}fprintf(somefile, "\n");}
		fclose(somefile);
		break;
		printf("I[K_AUG]...\t[DPOTRI_DRIVER]"
			"Attempt %d Status:%d, d %f.\n", i, info, d);
 	  if(info != 0 && try_d == 0){ /* Set-up stage */
			fprintf(stderr, "E[K_AUG]...\t[DPOTRI_DRIVER]"
				"info is non-zero ! %d\n", info);
			fprintf(stderr, "E[K_AUG]...\t[DPOTRI_DRIVER]"
				"Trying the modified Cholesky strategy\n");
		}
	  else if (info == 0){break;}
	  /* The previous fact was successful*/
	  else{dlast = d;}

  	if(dlast == 0.0 && try_d == 0){d = d0;}
  	else if(dlast == 0.0 && try_d > 0){d = kbp*d;}
  	else if(dlast > 0.0 && try_d > 0){d = kp*d;}
		else{d = (dmin > km*dlast) ? dmin:km*dlast;}

		if(d > dmax){
			fprintf(stderr, "E[K_AUG]...\t[DPOTRI_DRIVER]"
				"d > dmax ! %d\n", info);
			break;}
		for(j=0; j<n; j++){*(a + n * j + j) += d;} /* On all of them ?*/
		try_d++;
  }

	if(info != 0){
		fprintf(stderr, "W[K_AUG]...\t[DPOTRI_DRIVER]"
			"Maximum number of attempts reached ! %d\n", info);
		free(a);
		free(ai);
		free(aimax);
		exit(-1);
	}
	else{
		printf("I[K_AUG]...\t[DPOTRI_DRIVER]"
			"Cholesky fact, succesful!. %d\n", info);
	}
		/* Modified Cholesky heuristic
		Find the largest abs(a_ii) from the main diagonal
		Find the minimum & maximum off diagonal element
		Compute inf norm of column_i*/

	dpotri_(&uplo, &n, a, &lda, &info);


	if(info != 0){
		fprintf(stderr, "E[K_AUG]...\t[DPOTRI_DRIVER]"
		"info is non-zero ! %d", info);
		exit(-1);
	}
	else{
		printf("I[K_AUG]...\t[DPOTRI_DRIVER]"
		"Inversion , succesful!. %d\n", info);
	}

	somefile = fopen(_file_name_, "w");
	for(i=0; i<n; i++){
		for(j=0; j<n; j++){
			if(j>i){fprintf(somefile, "%.10g\t", *(a + n*j + i));}
			else{fprintf(somefile, "%.10g\t", *(a + n*i + j));}
			
		}
		fprintf(somefile, "\n");
	}
	fclose(somefile);
	free(a);
	free(ai);
	free(aimax);
}


double max_fun(int n, double *ai, int *imax){
	int i;
	double curr_max = 0;
	for(i=0; i<n; i++){
		if(fabs(ai[i]) > curr_max){
			curr_max = fabs(ai[i]);
			*imax = i;
		}
	}
	return curr_max;
}
