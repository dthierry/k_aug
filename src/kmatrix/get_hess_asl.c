#include <stdio.h>
#include <stdlib.h>
#include "get_hess_asl.h"


void get_hess_asl(ASL *asl, real *x, fint **Wcol, fint **Wrow, real **Wij, 
	int nvar, int ncon, int nobj, fint *n_nz_w, real *y, fint *nerror){

	real *c_body;
	real obj_val;
	real *Hcont; /* container for Hess */
	fint i, j, k;
	real ow;
	FILE *f_hess;
	/* setup the sparse hessian with multipliers */
	if (nobj == 0){
		ow = 0; /* set the objective weight to zero */
		*n_nz_w = sphsetup(0, ow, 1, 1);
		printf("I[KMATRIX]...\t[GET_HESS_ASL]"
			"No objective declared\n");
	}
	else{
		ow = 1; /* set the objective weight to one */
		*n_nz_w = sphsetup(-1, ow, 1, 1);
		printf("I[KMATRIX]...\t[GET_HESS_ASL]"
			"Objective found\n");
	}

	/* determine the kind of problem (min/max) */
	if (nobj > 0){
		if (nobj > 1){
			printf("W[KMATRIX]...\t[GET_HESS_ASL]"
			"The problem contains multiple obj_fun, will use the 1st one\n");
		}
		if(objtype[0]){
			printf("I[KMATRIX]...\t[GET_HESS_ASL]"
			"Maximization problem detected\n");
			/* set weight to -1 */
			ow = -1;}
	else{
		printf("I[KMATRIX]...\t[GET_HESS_ASL]"
			"Minimization problem detected\n");
		ow = 1;}
	}

	/* evaluate the objective function first  */
	obj_val = objval(0, x, nerror);

	c_body = (real *)malloc(sizeof(real) * ncon);
	/* need to evaluate the constraint-body first */
	conval(x, c_body, nerror);

	(*Wij) =  (real *)malloc(sizeof(real)*(*n_nz_w));
	(*Wcol) = (fint *)malloc(sizeof(fint)*(*n_nz_w));
	(*Wrow) = (fint *)malloc(sizeof(fint)*(*n_nz_w));

	Hcont = (real *)malloc(sizeof(real)*(*n_nz_w));

  f_hess = fopen("hess_debug.in", "w");

	/* Hessian of the Lagrange function matrix */
	if (*n_nz_w) {
	  if (n_obj > 0){
	    sphes(Hcont, -1, &ow, y);
	  }
	  else{
	    sphes(Hcont, 0, &ow, y);
  }


  /* pretty much compressed column format */
  k = 0;
  /* position the counter at 0 */
  for(i = 0; i < nvar; i++){
    for (j = sputinfo->hcolstarts[i]; j< sputinfo->hcolstarts[i+1]; j++){
      (*Wij)[k] = Hcont[j];
      (*Wcol)[k] = i + 1;
      (*Wrow)[k] = sputinfo->hrownos[j] + 1;
      fprintf(f_hess, "\t%ld\t%ld\t%.g\n", sputinfo->hrownos[j] + 1, i+1, Hcont[j]);
      k++;
    }
  }

	}
	fclose(f_hess);
	free(c_body);
	free(Hcont);
}