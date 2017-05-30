#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "get_hess_asl_aug.h"


void get_hess_asl_aug(ASL *asl, real *x, fint **Wcol, fint **Wrow, real **Wij, 
	int nvar, int ncon, int nobj, fint *n_nz_w, real *y, fint *nerror){
	int missing_nz;
	real *c_body;
	real obj_val;
	real *Hcont; /* container for Hess */
	fint i, j, k;
	real ow;
	FILE *f_hess;
	int *md_off_w;

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

	printf("I[KMATRIX]...\t[GET_HESS_ASL]"
	"Nonzeroes in the sparse hessian %d\n", (*n_nz_w));

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
	/* Added nvar to compensate for the case of not having elements in the main diagonal*/
	(*Wij) =  (real *)malloc(sizeof(real)*((*n_nz_w) + nvar));
	(*Wcol) = (fint *)malloc(sizeof(fint)*((*n_nz_w) + nvar));
	(*Wrow) = (fint *)malloc(sizeof(fint)*((*n_nz_w) + nvar));

	md_off_w = (int *)calloc(sizeof(int), nvar);

  f_hess = fopen("hess_debug.in", "w");

  k = 0;
  missing_nz = 0;
	/* Hessian of the Lagrange function matrix */
	if (*n_nz_w) {
		Hcont = (real *)malloc(sizeof(real)*(*n_nz_w));
	  if (n_obj > 0){
	    sphes(Hcont, -1, &ow, y);
	  }
	  else{
	    sphes(Hcont, 0, &ow, y);
  	}

  /* pretty much compressed column format */
  
  /* position the counter at 0 */
  for(i = 0; i < nvar; i++){
  	if(sputinfo->hcolstarts[i] - sputinfo->hcolstarts[i+1] == 0){
  		/* Missing column */
  		(*Wij)[k] = 0.0;
      (*Wcol)[k] = i + 1;
      (*Wrow)[k] = i + 1;
      fprintf(f_hess, "\t%ld\t%ld\t%.g\n", i+1, i+1, 0.0);
      md_off_w[i] = k++;
      missing_nz++;
  	}
  	else{
  		for (j = sputinfo->hcolstarts[i]; j< sputinfo->hcolstarts[i+1]; j++){
	      (*Wij)[k] = Hcont[j];
	      (*Wcol)[k] = i + 1;
	      (*Wrow)[k] = sputinfo->hrownos[j] + 1;
	      fprintf(f_hess, "\t%ld\t%ld\t%.g\n", sputinfo->hrownos[j] + 1, i+1, Hcont[j]);
	      k++;
	      }
    /* Last check */
      if((*Wcol)[k-1] != (*Wrow)[k-1]){
      	/* Not in the diagonal */
      	(*Wij)[k] = 0.0;
      	(*Wcol)[k] = i + 1;
      	(*Wrow)[k] = i + 1;
      	fprintf(f_hess, "\t%ld\t%ld\t%.g\n", i+1, i+1, 0.0);
      	md_off_w[i] = k++;
      	missing_nz++;
    	}
    	else{
    		md_off_w[i] = k-1;
    	}
  	}
  }
  free(Hcont);
	}
	else{
		printf("W[KMATRIX]...\t[GET_HESS_ASL]"
			"No non-zeroes in the hessian. Appending 0s at main-diag\n");
		/* No nz_ in the hessian; still want elements in md*/
		for (i = 0; i < nvar; i++)
		{
			(*Wij)[k] = 0.0;
			(*Wcol)[k] = i + 1;
			(*Wrow)[k] = i + 1;
			fprintf(f_hess, "\t%ld\t%ld\t%.g\n", i+1, i+1, 0.0);
			md_off_w[i] = k++;
		}
	}
	fclose(f_hess);
	free(c_body);
	printf("k-1 %d\n", k-1);
	printf("n_nz_w %d\n", (*n_nz_w));
	printf("n_nz_w+nvar %d\n", (*n_nz_w) + nvar);
	printf("missing_nz %d\n", missing_nz);
	assert(k-1 <= ((*n_nz_w) + nvar));
	f_hess = fopen("md_positions.in", "w");
	for(i=0; i< nvar; i++){
		fprintf(f_hess, "i %d\t%d\n", i, md_off_w[i]);
	}
	close(f_hess);
	
}