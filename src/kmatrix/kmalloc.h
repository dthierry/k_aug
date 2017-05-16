/* @source kmalloc.hfree(Wrow);
        free(Wcol);
        free(Wij);
**
** April 5th, 2017
** @author: David Thierry (dmolinat@andrew.cmu) dav0@lb2016-1

******************************************************************************

* @k_malloc ********************************************
**
** Read Calculate space for KKT matrix
**
** @param [r] Wrow
** @param [r] Wcol
** @return [fint] Space required
** @@
*******************************************************************************/
#ifndef KMALLOC_K
#define KMALLOC_K

#include "asl.h"
/* Calculates the amount of space for the whole KKT matrix */
fint k_malloc(fint *Wrow, fint *Wcol, real *Wij, fint Wnz, fint nvar, fint ncon,
 fint Anz); 
#endif /* KMALLOC_C */