
#include "ma57_driver.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "../../k_aug/inertia_strategy.h"


/*
   extern void ma57id_(double *CNTL, fint *ICNTL);
   extern void ma57ad_(fint *N, fint *NE, fint *IRN, fint *JCN, fint *LKEEP, fint *KEEP,
   fint *IWORK, fint *ICNTL, fint *INFO, double *RINFO);
   extern void ma57bd_(fint *N, fint *NE, double *A, double *FACT, fint *LFACT,
   fint *IFACT, fint *LIFACT, fint *LKEEP, fint *KEEP, fint *IWORK, fint *ICNTL,
   double *CNTL, fint *INFO, double *RINFO);
   extern void ma57cd_(fint *JOB, fint *N, double *FACT, fint *LFACT, fint *IFACT,
   fint *LIFACT,	fint *NRHS, double *RHS, fint *LRHS, double *WORK, fint *LWORK,
   fint *IWORK, fint *ICNTL, fint *INFO);
   extern void ma57dd_(fint *JOB, fint *N, fint *NE, double *A, fint *IRN, fint *JCN,
   double *FACT, fint *LFACT, fint *IFACT, fint *LIFACT, double *RHS, double *X,
   double *RESID, double *WORK, fint *IWORK, fint *ICNTL, double *CNTL,
   fint *INFO, double *RINFO);
   extern void ma57ed_(fint *N, fint *IC, fint *KEEP, double *FACT, fint *LFACT,
   double *NEWFAC, fint *LNEW, fint	*IFACT, fint *LIFACT, fint *NEWIFC, fint *LINEW,
   fint *INFO);*/

/*
   fint main(void){
   fint irn[] = {1, 2};
   fint jrn[] = {1, 2};
   double a[] = {1.0, 1.0};
   double b[] = {1.999999, 2.9e-1};
   double x[] = {0, 0};
   ma57_driver(2, 2, irn, jrn, a, 1, b, x);
   return 0;
   }

*/

/* MA57_DRIVER(fint *row_starts, fint *ia, fint *ja, double *a, fint n, fint n_rhs, double *b, double *x, fint nvar, fint ncon, fint no_inertia,
   fint nza, inertia_perts *inrt_pert, inertia_params inrt_parms, inertia_options *inrt_opts, double log10mu, linsol_opts ls_opts)
   */

void ma57_driver(fint *row_starts, fint *ia, fint *ja, double *a, fint n, fint n_rhs, double *b, double *x, fint nvar,
                 fint ncon, fint no_inertia,
                 fint nza, inertia_perts *inrt_pert, inertia_params inrt_parms, inertia_options *inrt_opts,
                 double log10mu, linsol_opts ls_opts) {

    double cntl[5];
    fint icntl[20];
    fint info[40];
    double rinfo[20];

    fint i, j;

    fint n_neig = 0;
    fint inertia_status = 1;
    fint reduce_pivtol;
    fint try_fact = 0;
    double trial_pivtol = ls_opts.pivot_tol0;
    double ratiorr = 0.0;
    fint inaccurateSol = 0;

    fint space_lk; /* for lkeep */
    fint *keep;
    fint *iwork;
    fint lfact = 0, lifact = 0;
    double *fact = NULL;
    fint *ifact = NULL;


    fint job;

    double *work;
    fint lwork;

    fint lrhs;


    double const residual_ratio_max = 1e-10;
    double *resid = NULL;

    fint incx = 1; /* for the norm calculation */
    double nrm_x = 0, nrm_r = 0;
    fint status = 0;
    printf("I[MA57]...\t[]"
           "***\n");

    space_lk = 5 * n + nza + (n > nza ? n : nza) + 42;
    lwork = n * n_rhs;

    keep = (fint *) malloc(sizeof(fint) * space_lk);
    iwork = (fint *) malloc(sizeof(fint) * 5 * n);
    resid = (double *) malloc(sizeof(double) * n);  /* the residuals */
    work = (double *) malloc(sizeof(double) * lwork);
    assert(keep);
    assert(iwork);
    assert(resid);
    assert(work);
    memset(keep, 0, sizeof(fint) * space_lk);
    memset(cntl, 0, sizeof(double) * 5);
    memset(icntl, 0, sizeof(fint) * 20);

    ma57id_(cntl, icntl);

    cntl[1 - 1] = 1e-08; /* initial pivtol*/
    icntl[5 - 1] = -1;
    icntl[11 - 1] = 16;
    icntl[12 - 1] = 16;
    ma57_analysis(&n, &nza, ia, ja, &space_lk, keep, iwork, icntl, info, rinfo);


    status = ma57_factorize(row_starts,
                            a, n,
                            nvar, ncon, no_inertia,
                            nza, inrt_pert, inrt_parms, inrt_opts,
                            log10mu, ls_opts,
                            &fact, &lfact, &ifact, &lifact, &space_lk,
                            keep, iwork, icntl, cntl, info, rinfo,
                            &reduce_pivtol, &trial_pivtol, &n_neig, &try_fact);

    if (status == 0) {
        printf("I[MA57]...\t[Factorize Success.]"
               "***\n");
    } else {
        printf("E[MA57]...\t[Factorize Failure!!!.]");
        exit(-1);
    }

    /**/
    /**/

    /* compute solution */
    ma57_solve(row_starts, a, ia, ja, n, nvar, ncon, nza, inrt_pert, inrt_parms, inrt_opts, log10mu, ls_opts, &fact,
               &lfact, &ifact, &lifact, &space_lk, keep, iwork, icntl, cntl, info, rinfo, work, &lwork, n_rhs, b, x,
               resid, &reduce_pivtol, &trial_pivtol, &n_neig, &try_fact, &ratiorr);

    n_neig = info[24 - 1];
    if (n_neig == ncon) {
        printf("W[K_AUG]...\t[MA57_DRIVER]"
               "Inertia check OK neig=%d, (neig == m).\n", n_neig);
    } else {
        printf("E[MA57]...\t[Inertia check Failure!!!.]");
        exit(-1);
    }
    /* actually i don't know if i have to reload b or x */

    free(keep);
    free(iwork);
    free(fact);
    free(ifact);
    free(work);
    free(resid);
}

fint ma57_analysis(fint *n, fint *nza, fint *ia, fint *ja, fint *space_lk, fint *keep, fint *iwork, fint *icntl, fint *info,
                  double *rinfo) {
    ma57ad_(n, nza, ia, ja, space_lk, keep, iwork, icntl, info, rinfo);
    if (info[0] != 0) {
        printf("I[MA57]...\t[ma57bd_]"
               "ma57ad_: info[0] is not zero; %d \n", info[1 - 1]);
        exit(-1);
    } else {
        printf("I[MA57]...\t[ma57bd_]"
               "ma57ad_: Status Ok.\n");
    }
    return 0;

}

fint ma57_factorize(const fint *row_starts, double *a, fint n, fint nvar, fint ncon, fint no_inertia, fint nza,
                   inertia_perts *inrt_pert, inertia_params inrt_parms, inertia_options *inrt_opts, double log10mu,
                   linsol_opts ls_opts, double **fact, fint *lfact, fint **ifact, fint *lifact, fint *space_lk, fint *keep,
                   fint *iwork, fint *icntl, double *cntl, fint *info, double *rinfo, fint *reduce_pivtol,
                   double *trial_pivtol, fint *n_neig, fint *try_fact) {
    fint i, j;
    fint lfact_new = *lfact, lifact_new = *lifact;
    fint ic = 0;
    fint inertia_status = 0;


    printf("I[MA57]...\t[MA57_FACTOR]\n");
    *lfact = (*lfact < info[9 - 1] * 2) ? info[9 - 1] * 2 : *lfact;
    *lifact = (*lifact < info[10 - 1] * 2) ? info[10 - 1] * 2 : *lifact;
    lfact_new = *lfact;
    lifact_new = *lifact;
    //if (*lifact < info[10 - 1] * 2){*lifact = info[10 - 1] * 2;}
    /* there is no way to know if fact and ifact have the appropriate size */
    *fact = (double *) (!*fact ? malloc(sizeof(double) * *lfact) : realloc(*fact, sizeof(double) * *lfact));
    *ifact = (fint *) (!*ifact ? malloc(sizeof(fint) * *lifact) : realloc(*ifact, sizeof(fint) * *lifact));
    assert(*fact);
    assert(*ifact);

    j = 0;
    if (inrt_opts->always_perturb_jacobian == 1) {
        printf("W[MA57]...\t[MA57_FACTOR]"
               "always_pert_jacobian is on before fact\n");
    }
    i = 0;
    while (i < ls_opts.max_inertia_steps) {
        /* factorization */
        if (info[0] == -3 || info[0] == -4) {
            /*MA57ED(N,IC,KEEP,FACT,LFACT,NEWFAC,LNEW,IFACT,LIFACT,NEWIFC,LINEW,INFO)
             * since realloc will copy the old arrays fact and ifact, we can re-use them.*/
            ma57ed_(&n, &ic, keep, *fact, lfact, *fact, &lfact_new, *ifact, lifact, *ifact, &lifact_new, info);
        } else {
            j = 0;
            ma57bd_(&n, &nza, a, *fact, lfact, *ifact, lifact, space_lk, keep, iwork, icntl, cntl, info, rinfo);
        }
        /* update these guys just in case */
        *lfact = lfact_new;
        *lifact = lifact_new;
        if (info[0] == -3 || info[0] == -4) { /* we need more memory */
            if (info[0] == -3) {
                ic = 0;
                lfact_new = info[17 - 1];
                *fact = (double *) realloc(*fact, sizeof(double) * lfact_new);
                assert(*fact);
            } else if (info[0] == -4) {
                ic = 1;
                lifact_new = info[18 - 1];
                *ifact = (fint *) realloc(*ifact, sizeof(fint) * lifact_new);
                assert(*ifact);
            } else {
                fprintf(stderr, "E[K_AUG]...\t[MA57_FACTOR]"
                                "Davs: We shouldn't reach this point!!\n");
                exit(-1);
            }
            /* check intertia? */
            i--;  /* This does not count  for the overall loop */
            j++;
            if (j > ls_opts.max_memory_al_steps) {
                fprintf(stderr, "E[K_AUG]...\t[MA57_FACTOR]"
                                "Reallocating Memory:Failed\n");
                exit(-1);
            }
            continue;  /* Try again. */
        } else { /* perhaps singular */
            if (info[0] == -10) { /* This one is problematic */
                fprintf(stderr, "E[K_AUG]...\t[MA57_DRIVER]"
                                "The KKT matrix is numerically singular. Assume delta_c > 0\n");
                fprintf(stderr, "E[K_AUG]...\t[MA57_DRIVER]"
                                "%d last known inertia\n", info[24 - 1]);
                if (inrt_pert->jacobian_perturbed == 1) {
                    fprintf(stderr, "E[K_AUG]...\t[MA57_DRIVER]"
                                    "Failure, the KKT matrix has been already perturbed\n");
                    exit(-1);
                }
                /*n_neig = id.info[12-1];*/
                *n_neig = 000000;
                inertia_status =
                        inertia_strategy(row_starts, a, nvar, ncon, *n_neig, inrt_pert, inrt_parms, inrt_opts, try_fact,
                                         log10mu,
                                         reduce_pivtol);
            }
            /*exit(-1);*/
        }
        /* check inertia */
        printf("I[K_AUG]...\t[MA57_FACTOR]"
               "n_neig = %d\n", info[24 - 1]);
        /* Get the number of negative eigenvalues */
        *n_neig = info[24 - 1];
        inertia_status =
                inertia_strategy(row_starts, a, nvar, ncon, *n_neig, inrt_pert, inrt_parms, inrt_opts, try_fact,
                                 log10mu,
                                 reduce_pivtol);

        if (inertia_status == 0) { return 0; }

        if (*reduce_pivtol != 0) {
            printf("I[K_AUG]...\t[MA57_FACTOR]"
                   "Increasing pivtol to fix inertia. pivot_tol %f\n", *trial_pivtol);
            cntl[1 - 1] = *trial_pivtol;
            *trial_pivtol = pow(*trial_pivtol, 0.75);
            if (*trial_pivtol >= ls_opts.pivtol_max) {
                fprintf(stderr, "E[K_AUG]...\t[MA57_FACTOR]"
                                "Failure, pivot tol is at it maximum.\n");
                fprintf(stderr, "W[K_AUG]...\t[MA57_FACTOR]"
                                "Inexact solution is now activated[Warning: results might not be good].\n");
                *trial_pivtol = ls_opts.pivtol_max;
                return 0;
                cntl[4 - 1] = 1e-08; /* static pivoting thingy */
            }
            /* Modify pivot tolerance */


        }
        i++;
    }


    return 1;
}


fint ma57_compute_ratios(const fint *row_start, const double *a, const fint *ja, fint n,
                        double *fact,
                        fint *lfact,
                        fint *ifact,
                        fint *lifact,
                        fint *iwork,
                        fint *icntl,
                        fint *info,
                        double *work,
                        fint *lwork,
                        fint n_rhs,
                        double *b,
                        double *x, double *resid,
                        double *ratiorr) {
    fint i, k;
    fint job = 1;
    fint incx = 1;
    double nrm_x = 1., nrm_r = 1.;
    for (i = 0; i < n * n_rhs; i++) {
        x[i] = b[i];
    }
    printf("I[MA57]...\t[MA57_RATIO]\n");
    ma57cd_(&job, &n, fact, lfact, ifact, lifact, &n_rhs, x, &n, work, lwork, iwork, icntl, info);
    nrm_x = dnrm2_(&n, x, &incx);

    if (info[0] != 0) {
        printf("W[MA57]...\t[ma57bd_]"
               "ma57bd_: info[0] is not zero; %d \n\n", info[1 - 1]);
        return 1;
    }

    memset(resid, 0, sizeof(double) * n);
    for (i = 0; i < n; i++) {
        for (k = row_start[i] - 1; k < row_start[i + 1] - 1; k++) {
            // printf("%d\t%d\t%d\n", i, k, ja[k]);
            resid[i] += a[k] * x[ja[k] - 1];
        }
        resid[i] -= b[i];
    }
    /* Residuals computation */
    nrm_r = dnrm2_(&n, resid, &incx);
    *ratiorr = nrm_r / (nrm_x + nrm_r);

    return 0;
}


fint ma57_solve(const fint *row_start, double *a, const fint *ia, const fint *ja, fint n, fint nvar,
               fint ncon, fint nza, inertia_perts *inrt_pert, inertia_params inrt_parms, inertia_options *inrt_opts,
               double log10mu, linsol_opts ls_opts,
               double **fact,
               fint *lfact,
               fint **ifact,
               fint *lifact, fint *space_lk, fint *keep,
               fint *iwork,
               fint *icntl, double *cntl,
               fint *info, double *rinfo,
               double *work,
               fint *lwork,
               fint n_rhs,
               double *b,
               double *x, double *resid,
               fint *reduce_pivtol, double *trial_pivtol, fint *n_neig, fint *try_fact, double *ratiorr) {

    fint inertia_status;
    fint inaccurateSol = 0;
    fint i, j;
    double ratio0; /* reference ratio */

    ma57_compute_ratios(row_start, a, ja, n, *fact, lfact, *ifact, lifact, iwork, icntl, info, work, lwork, n_rhs, b, x,
                        resid, ratiorr);
    printf("I[K_AUG]...\t[MA57_SOLVE]"
           ": Ratio of norm of scaled residuals (reported); %e \n", *ratiorr);
    ratio0 = *ratiorr;


    if (*ratiorr > ls_opts.residual_ratio_max) {
        printf("I[K_AUG]...\t[MA57_SOLVE]"
               ": The norm of residuals is larger than max ratio(computed)\n");
        inaccurateSol = 1;
    }

    if (inaccurateSol == 0) {
        printf("I[K_AUG]...\t[MA57_SOLVE]"
               "Accuracy at an acceptable level.\n\n");
        return 0;
    }
    /* Attempt to get a better solution. */
    j = 0;
    if (ls_opts.want_accurate == 1) {
        if (inaccurateSol == 1) {
            for (i = 0; i < ls_opts.max_refinement_steps; i++) {
                printf("I[K_AUG]...\t[MA57_SOLVE]"
                       ": Attempting to reduce residuals\n\n");
                if (inrt_pert->jacobian_perturbed == 0) {
                    printf("W[K_AUG]...\t[MA57_SOLVE]"
                           "Attempting to make dc > 0. (Jacobian regularization)\n");
                    inrt_opts->always_perturb_jacobian = 1;
                    inertia_status = inertia_strategy(ia, a, nvar, ncon, *n_neig, inrt_pert, inrt_parms, inrt_opts,
                                                      try_fact, log10mu,
                                                      reduce_pivtol); /*Perform correction*/
                } else {
                    printf("I[K_AUG]...\t[MA57_SOLVE]"
                           "Asking for better accuracy.\n");
                    *trial_pivtol = pow(*trial_pivtol, 0.5);
                    cntl[1 - 1] = *trial_pivtol;
                    /* Modify pivot tolerance */
                    if (*trial_pivtol > ls_opts.pivtol_max) {
                        fprintf(stderr, "E[K_AUG]...\t[MA57_SOLVE]"
                                        "Failure, pivot tol is at it maximum.\n");
                        break;
                    }

                }

                /*
                 * Analyze
                 * Factorize
                 *
                 * Try again if necessary */
                ma57_factorize(row_start, a, n, nvar, ncon, 0, nza, inrt_pert, inrt_parms, inrt_opts, log10mu, ls_opts,
                               fact, lfact, ifact, lifact, space_lk, keep, iwork, icntl, cntl, info, rinfo,
                               reduce_pivtol, trial_pivtol, n_neig, try_fact);
                if (info[0] < 0) {
                    printf("ERROR! STATUS RETURN: \tINFOG(1)= %d\n\t\t\t\tINFOG(2)= %d\n",
                           info[0], info[1]);
                    exit(-1);
                }
                /* compute solution */
                /*ratiorr = */
                ma57_compute_ratios(row_start, a, ja, n, *fact, lfact, *ifact, lifact, iwork, icntl, info, work, lwork,
                                    n_rhs, b, x, resid, ratiorr);
                printf("I[K_AUG]...\t[MA57_SOLVE]"
                       "Accuracy improvement tries %d, Trial pivtol %.g, Residual ratio %.g.\n", i, *trial_pivtol,
                       *ratiorr);
                if (*ratiorr < ls_opts.residual_ratio_max) {
                    printf("I[K_AUG]...\t[MA57_SOLVE]"
                           "Accuracy at acceptable level.\n\n");
                    break;
                } else if (fabs(ratio0 - *ratiorr) / ratio0 < 0.05) {
                    printf("I[K_AUG]...\t[MA57_SOLVE]"
                           "Accuracy is not improving.\n\n");
                    break;
                }
            }
        }
    }

    if (*ratiorr > ls_opts.residual_ratio_max) {
        fprintf(stderr, "W[K_AUG]...\t[MA57_DRIVER]\n\n"
                        "\t\tCould not fix the accuracy of the problem.\n"
                        "\t\tTry re-writing the problem or give a different point or change \"max_refinement_steps\"\n"
                        "\t\tWarning: results might be incorrect.\n"
                        "\t\tCurrent residual ratio %g; Max residual ratio %g.\n\n", *ratiorr,
                ls_opts.residual_ratio_max);
    }


}
