/* @source mumps_driver.c
**
** January 11th, 2018
** @author: David Thierry (dmolinat@andrew.cmu) dav0@lb2016-1

******************************************************************************

@pmumps_driver ********************************************
**
** interface to mumps, solves linear system
**
** @param [r] ia, sparse row
** @param [r] ja, sparse column
** @param [r] a, sparse matrix element
** @param [r] n, number of rows
** @param [r] nza, number of nz of a
** @param [r] nrhs, number of right hand sides
** @param [r] b, right hand side vector(matrix)
** @param [r] x, solution vector
** @param [r] nvar, number of variables
** @param [r] ncon, number of constraints
** @param [r] no_inertia, option to to deactivate the inertia check
** @param [r] nza, number of somethings
** @param [r] logmu0, logarithm of the barrier parameter.

** @@ This comes mostly from the example file at their website
*******************************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include "mumps_driver.h"
#include "mpi.h"
#include "../../../thirdparty/mumps512/MUMPS_5.1.2/include/dmumps_c.h"

/*int mumps_driver (fint *ia, fint *ja, real *a, fint n,
	fint n_rhs, real *b, real *x, fint nvar, fint ncon, int no_inertia, int nza, double logmu0);
 int mumps_driver(fint *ia, fint *ja, real *a, fi nt n,
	fint n_rhs, real *b, real *x, fint nvar, fint ncon, int no_inertia, int nza, double logmu0){*/
/*pardiso_driver(Kr_strt, Kcol, Kij, K_nrows, n_dof, rhs_baksolve, x_, n_var, n_con, no_inertia, nzK, logmu0, 1);*/
int mumps_driver(fint *ia, fint *ja, real *a, fint n,
	fint n_rhs, real *b, real *x, fint nvar, fint ncon, int no_inertia, int nza, double logmu0){
    /* Extra algorithmic steps:
     * inertia checking
     * d_r strategy
     * d_c strategy
     * soft-regularization
     * */
	DMUMPS_STRUC_C id;	/* datastructure */

    int myid, ierr;
    ierr = MPI_Init(NULL, NULL);
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	id.sym = 2;	/* symmetric indefinite matrix */
	id.comm_fortran = 2049;	/* i have no idea of what this is */
	id.par = 1;	/* the host is not involved */

	id.icntl[0] = -1;
	id.icntl[1] = -1;
	id.icntl[2] = -1;
	id.icntl[3] = 0;

	id.job = -1;  /* initial job */

    dmumps_c(&id);

    if (myid == 0){
        id.n = n;
        id.nnz = nza;
        id.irn = ia;
        id.jcn = ja;
        id.a = a;
        id.rhs = b;
    }

	id.n = n ;	/* number of rows */
	id.nnz = nza;	/* number of non-zeroes */
	id.irn = ia;	/* column coordinate */
	id.jcn = ja;	/* row coordinate */
	id.a = a;	/* matrix element */


	


	id.job = 1;  /* analysis */
    dmumps_c(&id);
	id.job = 2;	/* factorization */
    dmumps_c(&id);
	id.job = 3; /* compute solution */
    dmumps_c(&id);


	id.job = -2; /* terminate (deallocate) job */
    ierr = MPI_Finalize();
	return 0;
}