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
#ifndef ASSEM_RHSDCDP
#define ASSEM_RHSDCDP
#include "asl.h"

void assemble_rhs_dcdp(real **rhs_dcdp, fint nvar, fint ncon, int *n_p, int *n_x,
  SufDesc *dcdp, int **hr_point, SufDesc *var_order);
#endif /* ASSEM_RHSRH */