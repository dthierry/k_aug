/* @source untitled.c
** beta 01
** Month dayth, 20yy
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
#ifndef COMP_SIGMA
#define COMP_SIGMA

#include "../../thirdparty/asl/solvers/asl.h"
#define HUGE_NUMBER 1e300


void compute_sigma(ASL *asl, fint nvar, real *x, real *z_L, real *z_U, real *sigma, double logmu);

#endif /* COMP_SIGMA */