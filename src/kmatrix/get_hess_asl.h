#ifndef GET_HESS_ASL
#define GET_HESS_ASL
#include "asl.h"

void get_hess_asl(ASL *asl, real *x, fint **Wcol, fint **Wrow, real **Wij, 
	int nvar, int ncon, int nobj, fint *n_nz_w, real *y, fint *nerror);
#endif /*GET_HESS_ASL*/