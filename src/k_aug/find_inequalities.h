#ifndef FIND_INEQ_CON
#define FIND_INEQ_CON

#include "../../thirdparty/asl/solvers/asl.h"

void find_ineq_con(fint ncon_,real *LBC, int *c_flag);
void find_bounds  (fint nvar_,real *lbv);



#endif