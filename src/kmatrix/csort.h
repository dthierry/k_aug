#include <stdio.h>
#include <stdlib.h>
#include "asl.h"

typedef struct temp{
	fint *c;
	real *a;
}temp;

int compf(const void *t1, const void *t2);
int ccolum(fint *row, fint *col, real *aij, fint nz, fint nrow);

