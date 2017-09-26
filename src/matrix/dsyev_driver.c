/* @source untitled.c
** beta 01
** Month dayth, 20yy
** @author: David Thierry (dmolinat@andrew.cmu) dav0@lb2016-1

********************************************************************************

@fun_name ********************************************
**
** Driver for eigenvalue calculation
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

#include "dsyev_driver.h"
#include <stdio.h>
#include <stdlib.h>


int dsyev_driver(long n, double *_a, long Kn, int *sb_p){
	char jobz, uplo;
	double *w, *work;
	double *a;
	double w_mock[1], work_mock[1];
	int lwork, lda;
	int info=0;
	int i, j;
	int ret_val=0;
	FILE *somefile;
 /* transform a to up-triangular */
	
	jobz = 'N'; /* Eigenvalues only */
	uplo = 'U'; /* Uppertriangular */
	lda = n;
	lwork = -1;

	a = (double *)malloc(sizeof(double) * n * n);


	for(i=0; i<n; i++){
		for(j=0; j<n; j++){
			if(j > i){
				*(a + n * i + j) = 0.0;
			}
			else{
				*(a + n * i + j) = (*(_a + Kn * i + sb_p[j]) + *(_a + Kn * j + sb_p[i])) * 0.5;
			}
			
		}
	}

	dsyev_(&jobz, &uplo, &n, a, &lda, w_mock, work_mock,
	 &lwork, &info);

	printf("I[K_AUG]...\t[DSYEV_DRIVER]"
		"info %d, work %f\n", info, work_mock[0]);
	
	lwork = (int)work_mock[0];

	somefile = fopen("uptri.txt", "w");
	for(i=0; i<n; i++){
		for(j=0; j<n; j++){
			fprintf(somefile, "\t%.g", *(a + n * j + i));
		}
		fprintf(somefile, "\n");
	}
	fclose(somefile);

	w = (double * )calloc(sizeof(double), n);
	work = (double *)calloc(sizeof(double) , lwork);

	dsyev_(&jobz, &uplo, &n, a, &lda, w, work,
	 &lwork, &info);
	if(info != 0){
		printf("I[K_AUG]...\t[DSYEV_DRIVER]"
		"info is non-zero ! %d", info);
	}

	somefile = fopen("eig_red_hess.txt", "w");
	for(i=0; i<n; i++){
		fprintf(somefile, "\t%.g\n", w[i]);
	}
	fclose(somefile);

	/* Check Positive Definiteness*/
	for(i=0; i<n; i++){
		if(w[i] < 0.0){
			printf("I[K_AUG]...\t[DSYEV_DRIVER]"
				"Negative eigenvalue found i %d, val=%.6e.\n", i, w[i]);
		}
		if(w[i] < -1e-05){
			printf("I[K_AUG]...\t[DSYEV_DRIVER]"
				"The matrix is indefinite.");
			ret_val = 1;
			break;
		}
	}

	
	
	free(a);
	free(work);
	free(w);
	return ret_val;
}

