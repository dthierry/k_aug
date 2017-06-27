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


void dpotri_driver(int n, double *_a, long Kn, int *sb_p){
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
	dpotrf_(&uplo, &n, a, &lda, &info);

	if(info != 0){
		printf("I[KMATRIX]...\t[DPOTRI_DRIVER]"
		"info is non-zero ! %d", info);
		exit(-1);
	}
	else{
		printf("I[KMATRIX]...\t[DPOTRI_DRIVER]"
		"Cholesky fact, succesful!. %d\n", info);
	}

	/* dpotrf(uplo,n,a,lda,info) */
	/* dsyev_(&jobz, &uplo, &n, a, &lda, w, work,
	 &lwork, &info);*/
	dpotri_(&uplo, &n, a, &lda, &info);

	if(info != 0){
		printf("I[KMATRIX]...\t[DPOTRI_DRIVER]"
		"info is non-zero ! %d", info);
		exit(-1);
	}
	else{
		printf("I[KMATRIX]...\t[DPOTRI_DRIVER]"
		"Inversion , succesful!. %d\n", info);
	}

	somefile = fopen("inv_.txt", "w");
	for(i=0; i<n; i++){
		for(j=0; j<n; j++){
			if(j>i){fprintf(somefile, "%.10g\t", *(a + n*j + i));}
			else{fprintf(somefile, "%.10g\t", *(a + n*i + j));}
			
		}
		fprintf(somefile, "\n");
	}
	fclose(somefile);
	free(a);
}

