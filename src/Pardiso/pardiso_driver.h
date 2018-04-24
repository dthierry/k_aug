#ifndef PARDISO_DRIVER
#define PARDISO_DRIVER

#include "../ASL/solvers/asl.h"


extern void pardisoinit (void   *, int    *,   int *, int *, double *, int *);
extern void pardiso     (void   *, int    *,   int *, int *,    int *, int *, 
                  double *, int    *,    int *, int *,   int *, int *,
                     int *, double *, double *, int *, double *);
extern void pardiso_chkmatrix  (int *, int *, double *, int *, int *, int *);
extern void pardiso_chkvec     (int *, int *, double *, int *);
extern void pardiso_printstats (int *, int *, double *, int *, int *, int *,
                           double *, int *);

extern void pardiso_residual (int *mtype, int *n, double *a, int *ia, int *ja, double *b, double *x, double *y, double *normb, double *normr);

extern double dnrm2_(int *n, double *x, int *incx);

int pardiso_driver(fint *ia, fint *ja, real *a, fint n,
 fint nrhs, real *b, real *x, fint nvar, fint ncon, int no_inertia, int nza, double logmu0, char PardNeverGiveUp);

#endif /* PARDISO_DRIVER */
