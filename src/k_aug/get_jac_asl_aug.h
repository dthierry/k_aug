#ifndef GET_JAC_AUG
#define GET_JAC_AUG

#include "../../thirdparty/asl/solvers/asl.h"

typedef struct temp_v2{
	fint c;
	real a;
}temp_v2;

int compf_v2(const void *t1, const void *t2);
int get_jac_asl_aug(ASL *asl, real *x, fint *Acol, fint *Arow, real *Aij,
	int nvar, int ncon, fint nzc_, fint *nerror, int **nz_row_a);


#endif /* GET_JAC_AUG */