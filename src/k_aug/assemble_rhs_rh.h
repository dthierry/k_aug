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
#ifndef ASSEM_RHSRH
#define ASSEM_RHSRH
#include "../../thirdparty/asl/solvers/asl.h"

void assemble_rhs_rh(real **rhs_rh, fint nvar, fint ncon, int *n_dof, 
  SufDesc *var_f, int **hr_point);
#endif /* ASSEM_RHSRH */