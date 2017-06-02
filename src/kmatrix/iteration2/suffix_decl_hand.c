/* @source suffix_handler.c
** beta 01
** May 26th, 2017
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
#include "suffix_decl_hand.h"

#include <stdio.h>
#include <stdlib.h>


void suffix_decl_hand(int n_rhs, SufDecl *suf_ptr, char **ptr_name, char **rhs_name){
  static char sword[] = "rhs_";  
  int i;
  /*char dof_v[] = {"dof_v"};*/
  char second[] = {"con_flag"};

  /* First two suffixes */
  suf_ptr->name = ptr_name[0];
  suf_ptr->table = 0;
  suf_ptr->kind = ASL_Sufkind_var;
  suf_ptr->nextra = 0;

  (suf_ptr + 1)->name = ptr_name[1];
  (suf_ptr + 1)->table = 0;
  (suf_ptr + 1)->kind = ASL_Sufkind_var | ASL_Sufkind_outonly;
  (suf_ptr + 1)->nextra = 0;

  /* Suffix for the right hand side creation */
  

  /* this bit writes names of the rhs suffixes */
  if(n_rhs > 0){
    for(i=0; i < n_rhs; i++){
      char numeric_rhs[27];
      numeric_rhs[0] = '\0';
      sprintf(numeric_rhs, "%d", i); /* Number into string */
      *(rhs_name[i]) = '\0';
      strcat(rhs_name[i], sword);
      strcat(rhs_name[i], numeric_rhs);
    }

    for(i=0; i<n_rhs; i++){
      printf("NAME %s\n", rhs_name[i]);
    }

    for(i=0; i<n_rhs; i++){
      (suf_ptr + i + 2)->name = rhs_name[i];
      (suf_ptr + i + 2)->table = 0;
      (suf_ptr + i + 2)->kind = ASL_Sufkind_con|ASL_Sufkind_real;
      (suf_ptr + i + 2)->nextra = 0;
    }
  }

}
