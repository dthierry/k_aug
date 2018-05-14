//
// Created by dav0 on 5/2/18.
//
#include "config_kaug.h"
#include "mc19_driver.h"
#include <string.h>
#include <stdio.h>
//int mc30driver(fint n, fint nz, real *a, fint *irn, fint *icn, real *s)
int mc19driver(fint n, fint nz, real *a, fint *irn, fint *icn, real *s){
    real *w, *adummy, *r, *c;
    int i, nzdummy = 0;
    fint *irndummy, *icndummy;
    FILE *somefile;


    /* We know that it is not exactly 2 but it is good for now */
    adummy = (real *)malloc(2*nz * sizeof(real));
    irndummy = (fint *)malloc(2*nz * sizeof(fint));
    icndummy = (fint *)malloc(2*nz * sizeof(fint));
    r = (real *)malloc(n * sizeof(real));
    c = (real *)malloc(n * sizeof(real));

    for(i=0; i<nz; i++){
        if(irn[i] == icn[i]){
            adummy[nzdummy] = a[i];
            irndummy[nzdummy] = irn[i];
            icndummy[nzdummy] = icn[i];
            nzdummy++;
        }
        else{
            adummy[nzdummy] = a[i];
            irndummy[nzdummy] = irn[i];
            icndummy[nzdummy] = icn[i];
            nzdummy++;
            adummy[nzdummy] = a[i];
            irndummy[nzdummy] = icn[i];
            icndummy[nzdummy] = irn[i];
            nzdummy++;
        }
    }

    somefile = fopen("dummymatrix.txt", "w");
    for(i=0; i<nzdummy; i++){
        fprintf(somefile, "%d\t\t%d\t\t%.g\n", irndummy[i], icndummy[i], adummy[i]);
    }
    fclose(somefile);


     /* Working array */
    w = (real *)malloc(sizeof(real) * n * 5);
    memset(r, 0, n * sizeof(real));
    memset(c, 0, n * sizeof(real));

    /* Call mc19 */
    mc19ad_(&n, &nzdummy, adummy, irndummy, icndummy, r, c, w);

    /*Take the average and hope this works*/
    for(i = 0; i < n; i++){
        s[i] = (c[i] + r[i])/2;
    }

#ifndef PRINT_VERBOSE
    somefile = fopen("mc19_results.txt", "w");
    for(i=0; i<n; i++){
        fprintf(somefile, "%.g\t\t%.g\n", r[i], c[i]);
    }
    fclose(somefile);
#endif

    free(w);
    free(adummy);
    free(irndummy);
    free(icndummy);
    free(c);
    free(r);

    return 0;



}