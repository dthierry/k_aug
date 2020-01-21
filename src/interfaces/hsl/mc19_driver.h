//
// Created by dav0 on 5/2/18.
//

#ifndef K_AUG_MC19_DRIVER_H
#define K_AUG_MC19_DRIVER_H

#include "asl.h"

extern void mc19ad_(fint *nr, fint *nz, real *a, fint *irn, fint *icn, float *r, float *c, real *w);
int mc19driver(fint n, fint nz, real *a, fint *irn, fint *icn, real *r);

#endif /* K_AUG_MC19_DRIVER_H */

