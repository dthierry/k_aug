#ifndef GET_GRAD_F
#define GET_GRAD_F

#include "../../thirdparty/asl/solvers/asl.h"

void get_grad_f(ASL *asl, real *x, int nvar, real *gf, fint *nerror);

#endif /* GET_GRAD_F */