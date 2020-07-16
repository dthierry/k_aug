#ifndef MA57_DRIVER
#define MA57_DRIVER

#include "asl.h"
#include "../../k_aug/k_aug_data.h"

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
  fint *INFO);

extern double dnrm2_(fint *n, double *x, fint *incx);

void ma57_driver(fint *row_starts, fint *ia, fint *ja, double *a, fint n, fint n_rhs, double *b, double *x, fint nvar, fint ncon, fint no_inertia,
                 fint nza, inertia_perts *inrt_pert, inertia_params inrt_parms,
                 inertia_options *inrt_opts, double log10mu, linsol_opts ls_opts);

fint ma57_factorize(const fint *row_starts, double *a, fint n, fint nvar, fint ncon, fint no_inertia, fint nza,
                   inertia_perts *inrt_pert, inertia_params inrt_parms, inertia_options *inrt_opts, double log10mu,
                   linsol_opts ls_opts, double **fact, fint *lfact, fint **ifact, fint *lifact, fint *space_lk, fint *keep,
                   fint *iwork, fint *icntl, double *cntl, fint *info, double *rinfo, fint *reduce_pivtol,
                   double *trial_pivtol, fint *n_neig, fint *try_fact);

fint ma57_analysis(fint *n, fint *nza, fint *ia, fint *ja, fint *space_lk, fint *keep, fint *iwork, fint *icntl, fint *info,
                  double *rinfo);

fint ma57_compute_ratios(const fint *row_start, double const *a, const fint *ja, fint n,
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
                        double *x,
                        double *resid,
                        double *ratiorr);

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
               fint *reduce_pivtol, double *trial_pivtol, fint *n_neig, fint *try_fact, double *ratiorr);


#endif
