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
** todo: actually compute the norms and products with sparse linear algebra
** @@ This comes mostly from the example file at their website
*******************************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include "mumps_driver.h"

#include "../../../thirdparty/mumps/MUMPS/libseq/mpi.h"
#include "../../../thirdparty/mumps/MUMPS/include/dmumps_c.h"
#include "../../k_aug/inertia_strategy.h"

/*int mumps_driver (fint *ia, fint *ja, real *a, fint n,
	fint n_rhs, real *b, real *x, fint nvar, fint ncon, int no_inertia, int nza, double logmu0);
 int mumps_driver(fint *ia, fint *ja, real *a, fi nt n,
	fint n_rhs, real *b, real *x, fint nvar, fint ncon, int no_inertia, int nza, double logmu0){*/
/*pardiso_driver(Kr_strt, Kcol, Kij, K_nrows, n_dof, rhs_baksolve, x_, n_var, n_con, no_inertia, nzK, logmu0, 1);*/

/* we need the row-starts for inertia correction */


int
mumps_driver(fint *row_starts, fint *ia, fint *ja, double *a, fint n, int n_rhs, double *b, double *x, int nvar, int ncon, int no_inertia,
             int nza, inertia_perts *inrt_pert, inertia_params inrt_parms, inertia_options *inrt_opts, double log10mu, linsol_opts ls_opts) {
    /* Extra algorithmic steps:
     * inertia checking
     * d_r strategy
     * d_c strategy
     * soft-regularization
     * */
    DMUMPS_STRUC_C id;	/* datastructure */

    int myid, ierr;
    double *temp;
    int i, j;
    int n_neig=0;
    int inertia_status=1;
    int reduce_pivtol;
    int try_fact=0;
    double trial_pivtol = ls_opts.pivot_tol0;
    double ratiorr = 0.0;
    int inaccurateSol = 0;
    ierr = MPI_Init(NULL, NULL);
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    id.sym = 2;	/* symmetric indefinite matrix */
    id.comm_fortran = -654654;	/* i have no idea of what this is */
    id.par = 1;	/* the host is not involved */

    id.icntl[1-1] = -1;
    id.icntl[2-1] = -1;
    id.icntl[3-1] = -1;
    id.icntl[4-1] = -1;

    id.job = -1;  /* initial job */

    dmumps_c(&id);
    for(i=0; i< n*n_rhs; i++){x[i] = b[i];}
    if (myid == 0){
        id.n = n;
        id.nnz = nza;
        id.irn = ia;
        id.jcn = ja;
        id.a = a;
        id.rhs = b;
        id.nrhs = 1;
        id.lrhs = n;
        id.icntl[4 - 1] = -10;  /* Printing level */
        id.icntl[14 - 1] = 50; /* Memory factor */
        id.icntl[11 - 1] = 1;  /* Error analysis */
    }

    if(inrt_opts->always_perturb_jacobian == 1){fprintf(stderr, "always pert is on before fact\n");}

    /* inertia checking */
    j = 0;
    for(i=0; i<ls_opts.max_inertia_steps; i++){
        id.job = 1;  /* analysis */
        dmumps_c(&id);
        id.job = 2;	/* factorization */
        dmumps_c(&id);
        if(id.infog[0] == -9) {
            printf("I[K_AUG]...\t[MUMPS_DRIVER]"
                   "Reallocating Memory\n\n");
            id.icntl[14-1] = id.icntl[14-1] * 2 ;
            if(id.icntl[14-1] > 200){
                fprintf(stderr, "W[K_AUG]...\t[MUMPS_DRIVER]"
                                "icntl 14 > 200\n");
//                exit(-1);
            }
            i--;  /* This does not count  for the overall loop */
            j++;
            if (j > ls_opts.max_memory_al_steps) {
                fprintf(stderr, "E[K_AUG]...\t[MUMPS_DRIVER]"
                                "Reallocating Memory:Failed\n");
                exit(-1);
            }
            continue;  /* Try again. */
        }
        else if (id.infog[0]<0) {
            printf("ERROR! (PROC %d) STATUS RETURN: \tINFOG(1)= %d\n\t\t\t\tINFOG(2)= %d\n",
                   myid, id.infog[0], id.infog[1]);
            if (id.infog[0] == -10){
                ; /* This one is problematic */
            }
            exit(-1);
        }

        printf("I[K_AUG]...\t[MUMPS_DRIVER]"
               "n_neig = %d\n", id.info[12-1]);
        /* Get the number of negative eigenvalues */
        n_neig = id.info[12-1];
        inertia_status =
                inertia_strategy(row_starts, a, nvar, ncon, n_neig, inrt_pert, inrt_parms, inrt_opts, &try_fact, log10mu,
                                 &reduce_pivtol);

        if(inertia_status == 0){break;}

        if(reduce_pivtol != 0){
            printf("I[K_AUG]...\t[MUMPS_DRIVER]"
                   "Asking for better accuracy.\n");
            id.cntl[0] = trial_pivtol;
            /* Modify pivot tolerance */
            if(trial_pivtol > ls_opts.pivtol_max){
                fprintf(stderr, "E[K_AUG]...\t[MUMPS_DRIVER]"
                                "Failure, pivot tol is at it maximum.\n");
                exit(-1);
            }
            trial_pivtol = pow(trial_pivtol, 0.5);
        }

    }

    if(inertia_status == 1) {
        fprintf(stderr, "E[K_AUG]...\t[PARDISO_DRIVER]\n"
                        "\t\tCould not fix the inertia of the problem.\n"
                        "\t\tTry re-writing the problem or give a different point or change \"max_inertia_steps\"\n"
                        "\t\tError: exiting now.\n\n");
        exit(-1);
    }

    id.job = 3; /* compute solution */

    dmumps_c(&id);

    ratiorr = id.rinfog[6-1];
    /* Residuals computation */
    printf("I[K_AUG]...\t[MUMPS_DRIVER]"
           ": Ratio of norm of scaled residuals (reported); %e \n", ratiorr);

    if(ratiorr > ls_opts.residual_ratio_max){
        printf("I[K_AUG]...\t[MUMPS_DRIVER]"
               ": The norm of residuals is larger than max ratio(computed)\n");
        inaccurateSol = 1;
    }

    if(inaccurateSol == 0) {
        printf("I[K_AUG]...\t[MUMPS_DRIVER]"
               "Accuracy at an acceptable level.\n\n");
    }
    /* Attempt to get a better solution. */
    j = 0;
    if(ls_opts.want_accurate == 1){
        if(inaccurateSol == 1){
            for(i=0; i<ls_opts.max_refinement_steps; i++){

                printf("I[K_AUG]...\t[MUMPS_DRIVER]"
                       ": Attempting to reduce residuals\n\n");
                if(inrt_pert->jacobian_perturbed == 0) {
                    printf("W[K_AUG]...\t[MUMPS_DRIVER]"
                           "Attempting to make dc > 0. (Jacobian regularization)\n");
                    inrt_opts->always_perturb_jacobian = 1;
                    inertia_status = inertia_strategy(ia, a, nvar, ncon, n_neig, inrt_pert, inrt_parms, inrt_opts,
                                                      &try_fact, log10mu,
                                                      &reduce_pivtol); /*Perform correction*/
                }
                else{
                    printf("I[K_AUG]...\t[MUMPS_DRIVER]"
                           "Asking for better accuracy.\n");
                    id.cntl[0] = trial_pivtol;
                    /* Modify pivot tolerance */
                    if(trial_pivtol > ls_opts.pivtol_max){
                        fprintf(stderr, "E[K_AUG]...\t[MUMPS_DRIVER]"
                                        "Failure, pivot tol is at it maximum.\n");
                        exit(-1);
                    }
                    trial_pivtol = pow(trial_pivtol, 0.5);
                }

                id.job = 1;  /* analysis */
                dmumps_c(&id);
                id.job = 2;	/* factorization */
                dmumps_c(&id);
                if(id.infog[0] == -9) {
                    printf("I[K_AUG]...\t[MUMPS_DRIVER]"
                           "Reallocating Memory\n\n");
                    id.icntl[14-1] = id.icntl[14-1] * 2 ;
                    if(id.icntl[14-1] > 200){exit(-1);}
                    i--;  /* This does not count for the overall loop */
                    j++;
                    if (j > ls_opts.max_memory_al_steps) {
                        fprintf(stderr, "E[K_AUG]...\t[MUMPS_DRIVER]"
                                        "Reallocating Memory:Failed\n");
                    }
                    continue;  /* Try again. */
                }
                if (id.infog[0]<0) {
                    printf("ERROR! (PROC %d) STATUS RETURN: \tINFOG(1)= %d\n\t\t\t\tINFOG(2)= %d\n",
                           myid, id.infog[0], id.infog[1]);
                    exit(-1);
                }
                id.icntl[11 - 1] = 2;

                id.job = 3; /* compute solution */
                dmumps_c(&id);
                ratiorr = id.rinfog[6-1];

                printf("I[K_AUG]...\t[MUMPS_DRIVER]"
                       "Accuracy improvement tries %d, Perturbation val %.g, Resitual_ratio %.g.\n", i, trial_pivtol, ratiorr);
                if(ratiorr < ls_opts.residual_ratio_max){
                    printf("I[K_AUG]...\t[MUMPS_DRIVER]"
                           "Accuracy at acceptable level.\n\n");
                    break;
                }
            }
        }
    }

    if(ratiorr > ls_opts.residual_ratio_max){
        fprintf(stderr, "E[K_AUG]...\t[MUMPS_DRIVER]\n\n"
                        "\t\tCould not fix the accuracy of the problem.\n"
                        "\t\tTry re-writing the problem or give a different point or change \"max_refinement_steps\"\n"
                        "\t\tWarning: results might be incorrect.\n"
                        "\t\tCurrent residual ratio %g; Max residual ratio %g.\n\n", ratiorr, ls_opts.residual_ratio_max);
    }

    n_neig = id.info[12-1];
    if(n_neig==ncon){
        printf("W[K_AUG]...\t[MUMPS_DRIVER]"
               "Inertia check OK neig=%d, (neig == m).\n", n_neig);
    }
    id.nrhs = n_rhs;
    id.rhs = x;  /* Some weird behaviour happens if I don't do this */
    id.job = 3; /* compute solution */
    dmumps_c(&id);

    id.job = -2; /* terminate (deallocate) job */
    dmumps_c(&id);

    ierr = MPI_Finalize();

    return 0;
}
