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

#include "assemble_rhs_rh.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

void assemble_rhs_rh(real **rhs_rh, fint nvar, fint ncon, int *n_dof, 
  SufDesc *var_f, int **hr_point){
 	int j, i;
 	int temp;
 	FILE *somefile;
  int *help_list;

  if((nvar-ncon) <= 0){
    fprintf(stderr, "E[K_AUG]...\t[ASSM_RHS_RH]"
        "n_var <= n_con, exiting... \n");
    exit(-1);
  }

  help_list = (int*)malloc(sizeof(int)*(nvar-ncon));

 	
  /* Making the matrix by column as usual in fotran*/

 	(*n_dof) = 0; /* Count number of degrees of freedom */
  
  for(i=0; i < nvar; i++){
    temp = var_f->u.i[i]; /* Retrieve value of suffix*/
    /* Find non-zero */
    if(temp != 0){
      help_list[(*n_dof)] = i; /* Save the position of nz */
      (*n_dof)++;
      }
    }

  printf("I[K_AUG]...\t[ASSM_RHS_RH]"
    "According to the suffixes declared dof is %d \n", *(n_dof));
  if((*n_dof) > (nvar - ncon)){
    printf("E[K_AUG]...\t[ASSM_RHS_RH]"
      "No valid number of n_dof declared\n");
    exit(-1);
  }
  else if ((*n_dof) < (nvar - ncon) ){
    printf("I[K_AUG]...\t[ASSM_RHS_RH]"
      "n_dof is less than n_var - n_con \n");
  }
  else if ((*n_dof) == (nvar - ncon) ){
    printf("I[K_AUG]...\t[ASSM_RHS_RH]"
      "n_dof exactly n_var - n_con \n");
  }

  (*rhs_rh) = (real *)calloc((nvar + ncon) * (*n_dof), sizeof(real));
  (*hr_point) = (int *)malloc(sizeof(int)*(*n_dof));

  memset((*rhs_rh), 0, (nvar + ncon) * (*n_dof) * sizeof(real));

  for(i=0; i<(*n_dof); i++){
    (*hr_point)[i] = help_list[i];
    /*printf("i %d help %d\n", i, help_list[i]);*/
    (*rhs_rh)[(nvar + ncon) * i + (*hr_point)[i]] = 1.0;
  }
  free(help_list);

  somefile = fopen("rhs_red_hess", "w");
  
  
  for(i=0; i < (*n_dof); i++){
    fprintf(somefile, "\t%d\t", (*hr_point)[i]);
  }
  fprintf(somefile, "\n\n");

  for(i=0; i < (nvar+ncon); i++){
    for(j=0; j < (*n_dof); j++){
      fprintf(somefile, "\t%f\t", *((*rhs_rh) + j*(nvar+ncon) + i) );
    }
  fprintf(somefile, "\n");
  }


  fclose(somefile);

}
