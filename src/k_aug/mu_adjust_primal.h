/* @source untitled.c
** beta 01
** June 20th, 2017
** @author: David Thierry (dmolinat@andrew.cmu) dav0@lb2016-1

********************************************************************************

@fun_name ********************************************
**
** Computes mu and adjust the primal variables value if necesary.
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

#ifndef MU_ADJ
#define MU_ADJ

#include "../../thirdparty/asl/solvers/asl.h"

void mu_adjust_x(int nvar, double *x, double *lbv, real *zL, real *zU, double log10mu_target, double *logmu0);
#endif /*MU_ADJ*/
