#ifndef ASSEMBLE_RHSDSV3
#define ASSEMBLE_RHSDSV3

#include "asl.h"

void assemble_rhsds(int n_rhs, fint rhs_len, 
 real *rhsbksolv, fint nvar, fint ncon, SufDesc **rhs_ptr);

#endif /* ASSEMBLE_RHSDS */