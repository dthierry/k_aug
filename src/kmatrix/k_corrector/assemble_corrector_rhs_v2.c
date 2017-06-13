#include "assemble_corrector_rhs_v2.h"
#include "get_grad_f.h"
#include "get_jac_asl_aug.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

/* must traverse jac matrix by row */
/* @source k_assemble.c
**
** April 25th, 2017
** @author: David Thierry (dmolinat@andrew.cmu) dav0@lb2016-1

******************************************************************************

@assemble_corrector_rhs ********************************************
**
** Assembles [df + lamda*dC; cA] (dense)
**
** @param [r] Wrow
** @param [r] Wcol
** @param [r] Wij
** @@
*******************************************************************************/

void assemble_corrector_rhs_v2(ASL *asl, real *x, real *lambda,
  fint nvar, fint ncon,
	fint *Arow, fint *Acol, real *Aij, fint Anz,
  real **Crhs, fint *nerror, int *c_flag){
	int row, col;
	fint i;
	int j, l;
	/* column starts */

	/* gradient part */
	real *rhs;
	/* constraint part */
	real *rhs_con;
	/* full rhs */
	real *rhs_full;
	int *nz_row=NULL;

	/* A row and col transposes */
	
	FILE *somefile;
	memset(Arow, 0, sizeof(fint) * Anz);
	memset(Acol, 0, sizeof(fint) * Anz);
	memset(Aij,  0, sizeof(real) * Anz);

	get_jac_asl_aug(asl, x, Acol, Arow, Aij, nvar, ncon, Anz, nerror, &nz_row); /* Evaluate the gradient */
	if(*nerror != 0){
		fprintf(stderr, "Error in the evaluation of the Jacobian of the constraints \n");
	}
	assert(*nerror == 0);

	/* get gradient of f */
	rhs = (real *)malloc(sizeof(real) * (nvar));
	objgrd(0, x, rhs, nerror);
	if(*nerror != 0){
		fprintf(stderr, "Error in the evaluation of the Gradient of the objective. \n");
	}
	assert(*nerror == 0);
	/* get constraint body */
	rhs_con = (real *)malloc(sizeof(real) * (ncon));
	conval(x, rhs_con, nerror);
		if(*nerror != 0){
		fprintf(stderr, "Error in the evaluation of the constraint body. \n");
	}
	assert(*nerror == 0);

	/* traverse A */
	for(j=0; j<Anz; j++){
		row = Arow[j]-1;
		col = Acol[j]-1;
		assert(row < nvar);
		rhs[row] += lambda[col] * Aij[j];
	}

	rhs_full = (real *)calloc((ncon + nvar), sizeof(real));	

	for(i = 0; i<(nvar); i++){
		rhs_full[i] = -rhs[i];
	}

	j = 0;
	for(i = nvar; i<(nvar + ncon); i++){
		/* rhs_full[i] = rhs_con[j]; */
		rhs_full[i] = -rhs_con[j];
		if(c_flag[j] == -3){
			/* rhs_full[i] += -LUrhs[2*j]; original */
			rhs_full[i] += LUrhs[2*j];
		}
		else if(c_flag[j] == -1){
			/* rhs_full[i] += -LUrhs[2*j+1]; */
			rhs_full[i] += LUrhs[2*j+1];
		}
		else if(c_flag[j] == -2){
			/* rhs_full[i] += -LUrhs[2*j]; */
			rhs_full[i] += LUrhs[2*j];
		}
		
		j++;
	}
	somefile = fopen("rhs_corrector_unscaled.txt", "w");

	for(i=0; i<(nvar + ncon); i++){
		fprintf(somefile, "%d\t%f\n", i, rhs_full[i]);
	}
	fclose(somefile);

	free(nz_row);
	free(rhs);
	free(rhs_con);
	*Crhs = rhs_full;
}
