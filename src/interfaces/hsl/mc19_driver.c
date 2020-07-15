//
// Created by dav0 on 5/2/18.
//
#include "../../../config_kaug.h"
#include "mc19_driver.h"
#include <string.h>
#include <stdio.h>
#include <math.h>

/*int mc30driver(fint n, fint nz, real *a, fint *irn, fint *icn, real *s) */
int mc19driver(fint n, fint nz, real *a, fint *irn, fint *icn, real *s){
    double *w=NULL, *adummy=NULL;
    /* These guys must be single precision!! */
    float *r=NULL, *c=NULL;
    fint i, j=0, nzdummy = 0;
    fint *irndummy=NULL, *icndummy=NULL;
    FILE *somefile=NULL;

    /* We know that it is not exactly 2 but it is good for now */
    adummy = (double *)malloc(2*nz * sizeof(double));
    irndummy = (fint *)malloc(2*nz * sizeof(fint));
    icndummy = (fint *)malloc(2*nz * sizeof(fint));
    r = (float *)malloc(n * sizeof(float));
    c = (float *)malloc(n * sizeof(float));

    memset(adummy,   0, 2 * nz * sizeof(double));
    memset(irndummy, 0, 2 * nz * sizeof(fint));
    memset(icndummy, 0, 2 * nz * sizeof(fint));

    for(i=0; i<nz; i++){
        if(fabs(a[i]) < 1e-08){continue;}
        if(irn[i] == icn[i]){
            adummy[nzdummy] = a[i];
            irndummy[nzdummy] = irn[i];
            icndummy[nzdummy] = irn[i];
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

    if(nzdummy > 2*nz){
        printf("Bad memory allocation\n");
        exit(-1);
    }

    /*
    n = 4;
    nz = 8;
    nzdummy = 8;
    irndummy[0] = 4;
    irndummy[1] = 1;
    irndummy[2] = 4;
    irndummy[3] = 2;
    irndummy[4] = 3;
    irndummy[5] = 2;
    irndummy[6] = 3;
    irndummy[7] = 1;

    icndummy[0] = 4;
    icndummy[1] = 1;
    icndummy[2] = 2;
    icndummy[3] = 2;
    icndummy[4] = 1;
    icndummy[5] = 4;
    icndummy[6] = 3;
    icndummy[7] = 4;


    adummy[0] = 16000;
    adummy[1] = 100;
    adummy[2] = 14000;
    adummy[3] = 6;
    adummy[4] = 900;
    adummy[5] = 8;
    adummy[6] = 110000;
    adummy[7] = 4;*/



     /* Working array */
    w = (double *)malloc(sizeof(double) * n * 5);
    memset(r, 0, n * sizeof(float));
    memset(c, 0, n * sizeof(float));

    /* Call mc19 */
    mc19ad_(&n, &nzdummy, adummy, irndummy, icndummy, r, c, w);

    /*Take the average and hope this works*/
    for(i = 0; i < n; i++){
        s[i] = (c[i] + r[i])/2;
    }

#ifndef PRINT_VERBOSE
    somefile = fopen("mc19_results.txt", "w");
    for(i=0; i<n; i++){
        if(r[i] > 1e+4 || c[i] > 1e+4){
            printf("Warning: Scaling factor is too large.\n");
            j++;
        } 
        fprintf(somefile, "%.g\t\t%.g\n", r[i], c[i]);
    }
    fclose(somefile);
#endif

    if(j>0.01 * n){
        printf("Error MC19; set scaling factors to 0.\n");
        for(i = 0; i < n; i++){
            s[i] = 0.0;
        }
    }

    free(w);
    free(adummy);
    free(irndummy);
    free(icndummy);
    free(c);
    free(r);

    return 0;
}
