#ifndef WCREORD_H
#define WCREORD_H

#include "asl.h"

void k_assemble_cc(fint *Wrow, fint *Wcol, real *Wij, fint Wnz, fint nvar, fint ncon,
 fint *Arow, fint *Acol, real *Aij, fint Anz,
 fint *Krow, fint *Kcol, real *Kij, fint Knz,
 fint *Kr_strt);

#endif /* WCREORD_H */