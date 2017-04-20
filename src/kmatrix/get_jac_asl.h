#ifndef GET_JAC_ASL
#define GET_JAC_ASL

#include "asl.h"


void get_jac_asl(ASL *asl, real *x, fint *Acol, fint *Arow, real *Aij, 
	fint nzc_, fint *nerror);

#endif /* GET_JAC_ASL */