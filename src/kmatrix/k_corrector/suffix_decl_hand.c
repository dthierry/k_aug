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


void suffix_decl_hand(SufDecl *suf_ptr, char **ptr_name, char **rhs_name, unsigned n_reg_suf, int n_rhs){
  static char sword[] = "rhs_";  
  int i;
  
  /* First (Reg) suffixes */
  suf_ptr->name = ptr_name[0];
  suf_ptr->table = 0;
  suf_ptr->kind = ASL_Sufkind_var;
  suf_ptr->nextra = 0;

  (suf_ptr + 1)->name = ptr_name[1];
  (suf_ptr + 1)->table = 0;
  (suf_ptr + 1)->kind = ASL_Sufkind_var | ASL_Sufkind_outonly;
  (suf_ptr + 1)->nextra = 0;
  
  /* z_L bound mult */
  (suf_ptr + 2)->name = ptr_name[2];
  (suf_ptr + 2)->table = 0;
  (suf_ptr + 2)->kind = ASL_Sufkind_var | ASL_Sufkind_real;
  (suf_ptr + 2)->nextra = 0;
  /* z_U bound mult */
  (suf_ptr + 3)->name = ptr_name[3];
  (suf_ptr + 3)->table = 0;
  (suf_ptr + 3)->kind = ASL_Sufkind_var | ASL_Sufkind_real;
  (suf_ptr + 3)->nextra = 0;


  /* second_primal */
  (suf_ptr + 4)->name = ptr_name[4];
  (suf_ptr + 4)->table = 0;
  (suf_ptr + 4)->kind = ASL_Sufkind_var | ASL_Sufkind_real;
  (suf_ptr + 4)->nextra = 0;
  /* sedond_dual */
  (suf_ptr + 5)->name = ptr_name[5];
  (suf_ptr + 5)->table = 0;
  (suf_ptr + 5)->kind = ASL_Sufkind_con | ASL_Sufkind_real;
  (suf_ptr + 5)->nextra = 0;
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
      (suf_ptr + i + n_reg_suf)->name = rhs_name[i];
      (suf_ptr + i + n_reg_suf)->table = 0;
      (suf_ptr + i + n_reg_suf)->kind = ASL_Sufkind_con|ASL_Sufkind_real;
      (suf_ptr + i + n_reg_suf)->nextra = 0;
    }
  }

}
