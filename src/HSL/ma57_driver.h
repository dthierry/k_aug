#ifndef MA57_DRIVER
#define MA57_DRIVER

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

void ma57_driver(int n, int nz, int *i_rn, int *j_rn, double *a, int n_rhs,
	double *rhs, double *x);

#endif