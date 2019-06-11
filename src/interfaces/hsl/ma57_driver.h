#ifndef MA57_DRIVER
#define MA57_DRIVER

#include "../../../thirdparty/asl/solvers/asl.h"
#include "../../k_aug/k_aug_data.h"

extern void ma57id_(double *CNTL, int *ICNTL);
extern void ma57ad_(int *N, int *NE, int *IRN, int *JCN, int *LKEEP, int *KEEP,
 int *IWORK, int *ICNTL, int *INFO, double *RINFO);
extern void ma57bd_(int *N, int *NE, double *A, double *FACT, int *LFACT,
 int *IFACT, int *LIFACT, int *LKEEP, int *KEEP, int *IWORK, int *ICNTL,
 double *CNTL, int *INFO, double *RINFO);
extern void ma57cd_(int *JOB, int *N, double *FACT, int *LFACT, int *IFACT, 
	int *LIFACT,	int *NRHS, double *RHS, int *LRHS, double *WORK, int *LWORK,
	int *IWORK, int *ICNTL, int *INFO);
extern void ma57dd_(int *JOB, int *N, int *NE, double *A, int *IRN, int *JCN,
 double *FACT, int *LFACT, int *IFACT, int *LIFACT, double *RHS, double *X,
 double *RESID, double *WORK, int *IWORK, int *ICNTL, double *CNTL,
 int *INFO, double *RINFO);
extern void ma57ed_(int *N, int *IC, int *KEEP, double *FACT, int *LFACT,
 double *NEWFAC, int *LNEW, int	*IFACT, int *LIFACT, int *NEWIFC, int *LINEW,
  int *INFO);

extern double dnrm2_(int *n, double *x, int *incx);

void ma57_driver(fint *row_starts, fint *ia, fint *ja, double *a, fint n, int n_rhs, double *b, double *x, int nvar, int ncon, int no_inertia,
                 int nza, inertia_perts *inrt_pert, inertia_params inrt_parms,
                 inertia_options *inrt_opts, double log10mu, linsol_opts ls_opts);

int ma57_factorize(const fint *row_starts, double *a, fint n, int nvar, int ncon, int no_inertia, int nza,
                   inertia_perts *inrt_pert, inertia_params inrt_parms, inertia_options *inrt_opts, double log10mu,
                   linsol_opts ls_opts, double **fact, int *lfact, int **ifact, int *lifact, int *space_lk, int *keep,
                   int *iwork, int *icntl, double *cntl, int *info, double *rinfo, int *reduce_pivtol,
                   double *trial_pivtol, fint *n_neig, fint *try_fact);

int ma57_analysis(int *n, int *nza, int *ia, int *ja, int *space_lk, int *keep, int *iwork, int *icntl, int *info,
                  double *rinfo);

int ma57_compute_ratios(const fint *row_start, double const *a, const fint *ja, fint n,
                        double *fact,
                        int *lfact,
                        int *ifact,
                        int *lifact,
                        int *iwork,
                        int *icntl,
                        int *info,
                        double *work,
                        int *lwork,
                        fint n_rhs,
                        double *b,
                        double *x,
                        double *resid,
                        double *ratiorr);

int ma57_solve(const fint *row_start, double *a, const fint *ia, const fint *ja, fint n, int nvar,
               int ncon, int nza, inertia_perts *inrt_pert, inertia_params inrt_parms, inertia_options *inrt_opts,
               double log10mu, linsol_opts ls_opts,
               double **fact,
               int *lfact,
               int **ifact,
               int *lifact, int *space_lk, int *keep,
               int *iwork,
               int *icntl, double *cntl,
               int *info, double *rinfo,
               double *work,
               int *lwork,
               fint n_rhs,
               double *b,
               double *x, double *resid,
               int *reduce_pivtol, double *trial_pivtol, int *n_neig, fint *try_fact, double *ratiorr);


#endif
