/* @source mumps_driver.h
**
** January 11th, 2018
** @author: David Thierry (dmolinat@andrew.cmu) dav0@lb2016-1

******************************************************************************

@pmumps_driver ********************************************
**
** interface to mumps, solves linear system
**
** @param [r] ia, row starts
** @param [r] ja, sparse column
** @param [r] a, sparse matrix element
** @param [r] n, number of rows
** @param [r] nza, number of nz of a
** @param [r] nrhs, number of right hand sides
** @param [r] b, right hand side vector(matrix)
** @param [r] x, solution vector

*******************************************************************************/

#ifndef MUMPS_DRIVER
#define MUMPS_DRIVER

#include "asl.h"
#include "dmumps_c.h" /* header that contains types and mumps structure */

int mumps_driver (fint *ia, fint *ja, real *a, fint n, 
	fint n_rhs, real *b, real *x, fint nvar, fint ncon, int no_inertia, int nza, double logmu0);

#endif /* MUMPS_DRIVER */