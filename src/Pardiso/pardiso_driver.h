#ifndef PARDISO_DRIVER
#define PARDISO_DRIVER

#include "asl.h"


void pardisoinit (void   *, int    *,   int *, int *, double *, int *);
void pardiso     (void   *, int    *,   int *, int *,    int *, int *, 
                  double *, int    *,    int *, int *,   int *, int *,
                     int *, double *, double *, int *, double *);
void pardiso_chkmatrix  (int *, int *, double *, int *, int *, int *);
void pardiso_chkvec     (int *, int *, double *, int *);
void pardiso_printstats (int *, int *, double *, int *, int *, int *,
                           double *, int *);
int pardiso_driver(fint *ia, fint *ja, real *a, fint n, fint nza,
 fint nrhs, real *b, real *x);

#endif /* PARDISO_DRIVER */
