#ifndef ASS_CORR
#define ASS_CORR
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
#include "asl.h"
void assemble_corrector_rhs_v2(ASL *asl, real *x, real *lambda,
  fint nvar, fint ncon,
	fint *Arow, fint *Acol, real *Aij, fint Anz,
  real **Crhs, fint *nerror, int *c_flag);

#endif /* ASS_CORR */