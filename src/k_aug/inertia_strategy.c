/*
 * Created by dav0 on 4/24/18.
*/


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "inertia_strategy.h"
#include "k_aug_data.h"



int
inertia_strategy(int *row_strt, double *a, int nvar, int ncon, int n_eig, double *d_w, double *d_c, double *d_w_last,
                 double *d_c_last, inertia_params i_parm, inertia_options i_opts, int *try_n, double log10mu,
                 int *pert_pivot, int *jac_pert) {
    int j, k;
    double d_w_trial = 0.0;

    (*pert_pivot) = 0;
    if(i_opts.no_inertia){return 0;}

    if(n_eig > ncon){
        fprintf(stderr, "W[K_AUG]...\t[INERTIA_STRATEGY]"
                        "Wrong inertia(n_eig > m).\n");
        if((*d_w_last) == 0.0 && (*try_n) == 0){d_w_trial = i_parm.d_w0;}
        else if ((*d_w_last) == 0.0 && (*try_n) > 0){d_w_trial = i_parm.kbp * d_w_trial;}
        else if ((*d_w_last) > 0.0 && (*try_n) > 0){d_w_trial = i_parm.kp * d_w_trial;}
        else{d_w_trial = (i_parm.dmin > i_parm.km * (*d_w_last)) ? i_parm.dmin: i_parm.km * (*d_w_last);}

        if(d_w_trial > i_parm.dmax){
            fprintf(stderr, "E[K_AUG]...\t[INERTIA_STRATEGY]"
                            "The computed trial d_w is above maximum value.\n");
            exit(-1);
        }
        (*d_w) = d_w_trial - (*d_w);
        for(j=0; j<nvar; j++){
            k = row_strt[j]-1;
            a[k] += (*d_w); /* Add the difference */
        }
    }

    else if((n_eig < ncon)||(i_opts.always_perturb_jacobian)){
        fprintf(stderr, "W[K_AUG]...\t[PARDISO_DRIVER]"
                        "Wrong inertia(neig < m).\n");
        if((*jac_pert) == 0){
            (*jac_pert) = 1;
            fprintf(stderr, "W[K_AUG]...\t[PARDISO_DRIVER]"
                            "Attempting to make (*d_c) > 0.\n");
            (*d_c) = i_parm.dcb * pow((pow(10, log10mu)), i_parm.kc);
            for(j=nvar; j<nvar+ncon; j++){
                k = row_strt[j]-1;
                a[k] += -(*d_c);
            }
        }
        else{
            (*pert_pivot) = 1; /* Ask to change the pivot tolerance */
        }


        if((*jac_pert) == 1){
            fprintf(stderr, "W[K_AUG]...\t[PARDISO_DRIVER]"
                            "(*d_c) is already > 0\nThe Jacobian might be singular.\n");}
    }
        /*
         * else if(zi!=0){
         *         fprintf(stderr,"E[K_AUG]...\t[PARDISO_DRIVER]"
         *                                "Failure, there is a zero eigenvalue. The jacobian is possibly singular.\n");
         *                                        exit(-1);  Would this block ever get evaluated?}*/
    else{
        (*d_w_last) = d_w_trial;
        (*d_c_last) = (*d_c);
        printf("I[K_AUG]...\t[PARDISO_DRIVER]"
               "Inertia check successful.\n");
        return 0;
    }
    (*try_n)++;
    return 1;



}

