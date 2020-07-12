/* @source untitled.c
** beta 01
** Month dayth, 20yy
** @author: David Thierry (dmolinat@andrew.cmu) dav0@lb2016-1

********************************************************************************

@fun_name ********************************************
**
** Need better heuristics
** Description
**
** @param [r] 
** @param [r] 
** @param [r] 
** @param [r] 
** @param [r] 
** @param [r] 
** @param [r] 
** @return something
*******************************************************************************/

#include "sigma_compute.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "../../config_kaug.h"

#define HUGE_NUMBER 1e300


void compute_sigma(ASL *asl, fint nvar, real *x, real *z_L, real *z_U, real *sigma, double logmu){
    int i, j=0;
    real sl=0, su=0;
    double machine_epsi = DBL_EPSILON;
    /*FILE *s_file; */
    /*
    if(!(suf_zL->u.r)){
        fprintf(stderr, "W[K_AUG]...\t[K_AUG_ASL]"
        "No ipopt_zL_out suffix declared, setting zL = 0.\n");

    }
    else{
        for(i=0; i< nvar; i++){
            z_L[i] = suf_zL->u.r[i];
        }
    }
    if(!(suf_zU->u.r)){
        fprintf(stderr, "W[K_AUG]...\t[K_AUG_ASL]"
        "No ipopt_zU_out suffix declared, setting zU = 0.\n");

    }
    else{
        for(i=0; i< nvar; i++){
            z_U[i] = suf_zU->u.r[i];
        }
    }
    */
    for(i=0; i<nvar; i++){
        if(LUv[2*i] < -HUGE_NUMBER){
            sl = 0.0;
        }
        else if (fabs(x[i] - LUv[2 * i]) < machine_epsi * exp(logmu)){
            fprintf(stderr, "E[K_AUG]...\t[SIGMA_COMPUTE]"
                            "is x[%d] = xlb?\n", i);
            j++; /* Should we fix sl or su to something else instead:? */
        }
        else{
            sl = z_L[i]/(x[i] - LUv[2*i]);
        }
        /**/
        if(LUv[2*i + 1] > HUGE_NUMBER){
            su = 0.0;
        }
        else if (fabs(-x[i] + LUv[2 * i + 1]) < machine_epsi * exp(logmu)){
            fprintf(stderr, "E[K_AUG]...\t[SIGMA_COMPUTE]"
                            "is x[%d] = xlb?\n", i);
            j++; /* Should we fix sl or su to something else instead:? */
        }
        else{
            /* now with correct expression */
            su = z_U[i]/(x[i] - LUv[2*i+1]);
        }
        sigma[i] = sl + su;
    }

#ifndef PRINT_VERBOSE
    s_file = fopen("sigma_out", "w");
    fprintf(s_file, "%s\t%s\t%s\t%s\t%s\t%s\t\t%s\n",
            "i", "z_L[i]", "z_U[i]", "LUv[2*i]", "x[i]", "LUv[2*i+1]", "sigma[i]");
    for(i=0; i<nvar; i++){
        fprintf(s_file, "%d\t%.g\t%.g\t%.g\t%.g\t%.g\t\t%e %c\n",
                i, z_L[i], z_U[i], LUv[2*i], x[i], LUv[2*i+1], sigma[i], sigma[i]>1e-02 ? '*':' ');
    }
    fclose(s_file);
#endif
    if(j > 0){
        fprintf(stderr, "E[K_AUG]...\t[SIGMA_COMPUTE]"
                        "Unresolved sigma values, is the mu computation successful? or is this a square problem?\n"
                        "Check the sigma_out file with the PRINT_VERBOSE (cmake) option set to 0\n");
        /*exit(-1);*/
    }
}

