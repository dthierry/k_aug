#ifndef WNZAPP_H
#define WNZAPP_H
#include "asl.h"
void wnzappnd(fint *Wrow, fint *Wcol, real *Wij, fint Wnz, fint nvar,
	fint *Wr_new, fint *Wc_new, real *Wi_new, fint Wnz_new);
#endif /* WNZAPP_H */