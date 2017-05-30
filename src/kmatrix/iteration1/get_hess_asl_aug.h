#ifndef HESS_AUG
#define HESS_AUG

#include "asl.h"


void get_hess_asl_aug(ASL *asl, real *x, fint **Wcol, fint **Wrow, real **Wij, 
	int nvar, int ncon, int nobj, fint *n_nz_w, real *y, fint *nerror);

#endif /* HESS_AUG */