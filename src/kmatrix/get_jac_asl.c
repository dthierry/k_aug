#include <stdio.h>
#include <stdlib.h>
#include "get_jac_asl.h"


void get_jac_asl(ASL *asl, real *x, fint *Acol, fint *Arow, real *Aij,
 fint nzc_, fint *nerror){
	int j;
	cgrad *cg, **cgx;
	FILE *f_jac;
	real Jcont[nzc_]; // container for Jac
	int i;
	j = 0;
	// Jacobian
	// check if there are non-zeroes in the Jacobian
	if(nzc_ <= 0){
		printf("[KMATRIX]...\t[GET_JAC_ASL]"
			"The jacobian has no structural non-zeroes; exiting\n");
		exit(1);
	}
	
	// evaluates the Jacobian matrix
	jacval(x, Jcont, nerror);

	if(*nerror != 0){
		printf("[KMATRIX]...\t[GET_JAC_ASL]nerror points to %ld\n", *nerror);
		exit(1);
	}
	f_jac = fopen("jacobi_debug.in", "w");

	cgx = Cgrad; // point to the beggining of the list

	for(i = 1; i <= n_con; i++) {
    // moves by constraint
    if((cg = *cgx++)) {
    // iterates for a given constraint
	    do{
	    // moves by nz in the constraint
		    fprintf(f_jac, "%d\t%d\t%.g\n",i , cg->varno+1, Jcont[cg->goff]);
	      Arow[j] = i;
	      Acol[j] = cg->varno+1;
	      Aij[j] = Jcont[cg->goff];
	      j++;
		    }
	    while ((cg = cg->next)) ;
	  }
  }

	fclose(f_jac);

}