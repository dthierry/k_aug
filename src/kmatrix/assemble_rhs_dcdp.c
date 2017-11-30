/* @source assemble_rhsds_red_hess.c
** beta 0
** April 18th, 2017
** @author: David Thierry (dmolinat@andrew.cmu) dav0@lb2016-1

********************************************************************************

@assemble_rhsds_red_hess ********************************************
**
** Assembles right hand sides for the reduced hessian
**
** @param [r] n_rhs
** @param [r] rhs_len
** @param [r] rhsbksolv
** @param [r] dp
** @param [r] nvar
** @param [r] ncon
** @param [r] rhs_ptr
** @@
*******************************************************************************/

#include "assemble_rhs_dcdp.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

void assemble_rhs_dcdp(real **rhs_dcdp, fint nvar, fint ncon, int *n_p, 
  SufDesc *dcdp, int **hr_point){
 	int j, i;
 	int temp;
 	FILE *somefile;
  int *help_list;

  /* We don't know the amount of rhs's so allocate the maximum possible namely n+m*/
  help_list = (int*)malloc(sizeof(int)*(nvar+ncon));
  /* Then use hr_point to mark the positions (KKT) indexed by the corresponding rhs*/
 	
  /* Making the matrix by column as usual in fotran*/

 	(*n_p) = 0; /* Count number of parameters */
  
  for(i=0; i < ncon; i++){
    temp = dcdp->u.i[i]; /* Retrieve value of suffix*/
    /* Find non-zero */
    if(temp != 0){
      help_list[(*n_p)] = i + nvar; /* Save the position of nz */
      (*n_p)++;
      }
    }

  printf("I[KMATRIX]...\t[ASSM_RHS_DCDP]"
    "According to the suffixes declared len p is %d \n", *(n_p));
  if((*n_p) <= (0)){
    printf("E[KMATRIX]...\t[ASSM_RHS_DCDP]"
      "No valid number of n_p declared\n");
    exit(-1);
  }
  

  (*rhs_dcdp) = (real *)calloc((nvar + ncon) * (*n_p), sizeof(real));
  (*hr_point) = (int *)malloc(sizeof(int)*(*n_p));

  memset((*rhs_dcdp), 0, (nvar + ncon) * (*n_p) * sizeof(real));

  for(i=0; i<(*n_p); i++){
    (*hr_point)[i] = help_list[i];
    /*printf("i %d help %d\n", i, help_list[i]);*/
    (*rhs_dcdp)[(nvar + ncon) * i + (*hr_point)[i]] = 1.0;
  }
  free(help_list);

  somefile = fopen("rhs_dcdp", "w");
  
  
  for(i=0; i < (*n_p); i++){
    fprintf(somefile, "\t%d\t", (*hr_point)[i]);
  }
  fprintf(somefile, "\n\n");

  for(i=0; i < (nvar+ncon); i++){
    for(j=0; j < (*n_p); j++){
      fprintf(somefile, "\t%f\t", *((*rhs_dcdp) + j*(nvar+ncon) + i) );
    }
  fprintf(somefile, "\n");
  }


  fclose(somefile);

}
