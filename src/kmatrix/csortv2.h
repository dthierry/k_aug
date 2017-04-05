#ifndef CCSORT_H
#define CCSORT_H
#include "asl.h"

typedef struct temp{
	fint *c;
	real *a;
}temp;
	
int compf(const void *t1, const void *t2);
int sortcol(fint *row, fint *col, real *aij, fint nz, fint nrow);

#endif /* CCSORT_H */