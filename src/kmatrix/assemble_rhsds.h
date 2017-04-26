#ifndef ASSEMBLE_RHSDS
#define ASSEMBLE_RHSDS

#include "asl.h"

void assemble_rhsds(int n_rhs, fint rhs_len, 
 real **rhsbksolv, fint nvar, fint ncon, SufDesc *rhs_ptr);

#endif /* ASSEMBLE_RHSDS */