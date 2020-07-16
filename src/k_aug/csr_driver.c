/* @source assemble_kkt_by_row.c
** beta 01
** May 31st, 2017
** @author: David Thierry (dmolinat@andrew.cmu) dav0@lb2016-1

********************************************************************************

@fun_name ********************************************
**
** Stores the KKT matrix in the following way:
** Similar to CSR, we have arrays Kcol and Kij for the non-zeroes on each row
** Kr_strt contains the rowstarts for a given row
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

#include "csr_driver.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "../../config_kaug.h"

void csr_driver(int nvar, int ncon, int nzW, int nzA,
	int *nzr_w, int *nzr_a,
	fint *Wrow, fint *Wcol, real *Wij,
	fint *Arow, fint *Acol, real *Aij,
	fint **Kr, fint **Kc, real **K, fint **Kr_strt){
	
	int i, j, k;
	int row;
	int *nzr_k;
	int **Kcol;
	real **Kij;
	int *rn;
	int sum_rn=0;
	FILE *somefile;
	int Pardiso_flag=1; /* For the moment */

	fint nzK;


	nzr_k = (int *)malloc(sizeof(int) *(nvar+ncon));
	rn =    (int *)calloc(sizeof(int), (nvar+ncon)); /* This better be 0's*/
	
	assert(nzr_k != NULL);
	assert(rn != NULL);

	for(i=0; i < nvar; i++){
		nzr_k[i] = nzr_w[i] + nzr_a[i]; /* nz per row */
		/*fprintf(stderr, "i %d nzrk %d \n", i, nzr_k[i]);*/
	}

	Kcol = (int **)malloc(sizeof(int *) * (nvar+ncon));
	Kij  = (real**)malloc(sizeof(real *) * (nvar+ncon));
	assert(Kcol != NULL);
	assert(Kij  != NULL);
	
	*Kr_strt = (fint *)malloc(sizeof(int *) * (nvar+ncon+1));
	assert(*Kr_strt != NULL);
	
	
	for(i=0; i < nvar; i++){
		Kcol[i] = (int *)malloc(sizeof(int) * nzr_k[i]);
		Kij[i]  = (real*)malloc(sizeof(real) *nzr_k[i]);
	}

	
	/* Traverse through the hessian of the lagrangian*/
	for(i=0; i < nzW; i++){
		row = (Wrow[i]-1);
		assert(rn[row] < nzr_k[row]);

		Kcol[row][rn[row]] = Wcol[i];
		Kij [row][rn[row]] = Wij[i];
				
		rn[row]++;
	}

	/* Traverse through the Gradients of the constraints*/
	for(i=0; i < nzA; i++){
		row = Arow[i]-1;
		assert(rn[row] <= nzr_k[row]);

		Kcol[row][rn[row]] = Acol[i] + nvar;
		Kij [row][rn[row]] = Aij[i];
		
		rn[row]++;
	}

	/* Note that the first element of each K array is on the main
	diagonal */

	/* If Pardiso is the linear solver
	 * Trailing zeroes of the lower bottom right of the kkt matrix*/
	if(Pardiso_flag){
		for(i=nvar; i < (ncon + nvar); i++){
			Kcol[i] = (int *)malloc(sizeof(int));
			Kij[i]  = (real*)malloc(sizeof(real));
		}
		for(i=nvar; i < (nvar+ncon); i++)	{
			Kcol[i][0] = i + 1;
			Kij [i][0] = 0E+00;
			rn[i]++;	
		}	
	}
#ifndef PRINT_VERBOSE
	somefile = fopen("row_start.in", "w");
#endif
	(*Kr_strt)[0] = 1;
#ifndef PRINT_VERBOSE
	fprintf(somefile, "\t%d\n", (*Kr_strt)[0]);
#endif
	sum_rn += rn[0] + 1;
	
	for(i=1; i<nvar; i++){
		(*Kr_strt)[i] = sum_rn;
#ifndef PRINT_VERBOSE
		fprintf(somefile, "\t%d\n", (*Kr_strt)[i]);
#endif
		sum_rn += rn[i];
	}

	if(Pardiso_flag){
		for(i=nvar; i < (ncon + nvar); i++){
			(*Kr_strt)[i] = sum_rn;
#ifndef PRINT_VERBOSE
			fprintf(somefile, "\t%d\n", (*Kr_strt)[i]);
#endif
			sum_rn += rn[i];
		}
		(*Kr_strt)[(ncon + nvar)] = sum_rn;
#ifndef PRINT_VERBOSE
		fprintf(somefile, "\t%d\n", (*Kr_strt)[(ncon + nvar)]);
#endif
	}
	/*else do something when not pardiso*/

#ifndef PRINT_VERBOSE
	fclose(somefile);
#endif
	/* If lower triangular*/
	if(Pardiso_flag){nzK = nzA + nzW + ncon;}
	else{nzK = nzA + nzW;}

	*Kr = (fint *)malloc(sizeof(fint)*nzK);
	*Kc = (fint *)malloc(sizeof(fint)*nzK);
	*K  = (real *)malloc(sizeof(real)*nzK);

	assert(*Kr != NULL);
	assert(*Kc != NULL);
	assert(*K  != NULL);

	k = 0;
#ifndef PRINT_VERBOSE
	somefile = fopen("kkt.in", "w");
#endif
	for(i=0; i<(nvar + ncon); i++){
		for(j=0; j<rn[i]; j++){
#ifndef PRINT_VERBOSE
            fprintf(somefile, "\t%d\t%d\t%.g\n", i+1, Kcol[i][j], Kij[i][j]);
#endif
		(*Kr)[k] = i + 1;
		(*Kc)[k] = Kcol[i][j];
		(*K)[k]  = Kij[i][j];
		k++;
		}
	}
#ifndef PRINT_VERBOSE
	fclose(somefile);
#endif
	
	assert(k == nzK);

	for(i=0; i<nvar; i++){
		assert(i+1 == Kcol[i][0]);
	}


	for(i=0; i < (nvar+ncon); i++){
		free(Kcol[i]);
		free(Kij[i]);
	}
	free(Kcol);
	free(Kij);
	free(rn);
	free(nzr_k);
}

