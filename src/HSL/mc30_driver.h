/* @source mc30_driver.c
**
** April 10th, 2017
** @author: David Thierry (dmolinat@andrew.cmu) dav0@lb2016-1

******************************************************************************

@mc30driver ********************************************
**
** Calculates scaling factors
**
** @param [r] nr row number
** @param [r] nz number of nz
** @@
*******************************************************************************/
#ifndef _MC30DRIVER_
#define _MC30DRIVER_


#include "../../thirdparty/asl/solvers/asl.h"

extern void mc30ad_(fint *nr, fint *nz, real *a, fint *irn, fint *icn, real *s, real *w, fint *lp, fint *ifail);

int mc30driver(fint n, fint nz, real *a, fint *irn, fint *icn, real *s);

#endif /* _MC30DRIVER */