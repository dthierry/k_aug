/*
 * Created by dav0 on 4/24/18.
*/
/*
 * @param row_strt*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "inertia_strategy.h"
#include "k_aug_data.h"




int
inertia_strategy(int *row_strt, double *a, int nvar, int ncon, int n_eig, inertia_perts *i_pert, inertia_params i_parm,
                 inertia_options *i_opts, int *try_n, double log10mu, int *pert_pivot) {
    int j, k;
    double d_w_trial = i_pert->d_w;

    (*pert_pivot) = 0; /* Initialize this little guy always. */
    if(i_opts->no_inertia){return 0;}

    /*
     *
     */
    /*
    printf("within inertiastrat kp %g\n", i_parm.kp);
    printf("within inertiastrat kbp %g\n", i_parm.kbp);
    printf("within inertiastrat kc %g\n", i_parm.kc);
    printf("within inertiastrat dcp %g\n", i_parm.dcb);
    printf("within inertiastrat d_w0 %g\n", i_parm.d_w0);
    printf("within inertiastrat dmin %g\n", i_parm.dmin);
    printf("within inertiastrat dmax %g\n", i_parm.dmax);
    */
    /* always check this first */
    if(n_eig > ncon){
        fprintf(stderr, "W[K_AUG]...\t[INERTIA_STRATEGY]"
                        "Wrong inertia(n_eig > m).\n");
        if(i_pert->d_w_last == 0.0 && (*try_n) == 0){d_w_trial = i_parm.d_w0;}
        else if (i_pert->d_w_last == 0.0 && (*try_n) > 0){d_w_trial = i_parm.kbp * d_w_trial;}
        else if (i_pert->d_w_last > 0.0 && (*try_n) > 0){d_w_trial = i_parm.kp * d_w_trial;}
        else{d_w_trial = (i_parm.dmin > i_parm.km * i_pert->d_w_last) ? i_parm.dmin: i_parm.km * i_pert->d_w_last;}

        if(d_w_trial > i_parm.dmax){
            fprintf(stderr, "E[K_AUG]...\t[INERTIA_STRATEGY]"
                            "The computed trial d_w is above maximum value.\n");
            exit(-1);
        }
        i_pert->d_w = d_w_trial - i_pert->d_w;

        for(j=0; j<nvar; j++){
            k = row_strt[j]-1;
            a[k] += i_pert->d_w; /* Add the difference */
        }
    }

    else if((n_eig < ncon)||(i_opts->always_perturb_jacobian == 1)){
        if(i_opts->always_perturb_jacobian == 1){
            fprintf(stderr, "W[K_AUG]...\t[INERTIA_STRATEGY]"
                            "Always perturb jacobian is on.\n");
        }
        else{fprintf(stderr, "W[K_AUG]...\t[INERTIA_STRATEGY]"
                             "Wrong inertia(neig < m).\n");
        }


        if(i_pert->jacobian_perturbed == 0){
            i_pert->jacobian_perturbed = 1;
            fprintf(stderr, "W[K_AUG]...\t[INERTIA_STRATEGY]"
                            "Attempting to make i_pert->d_c > 0.\n");
            i_pert->d_c = i_parm.dcb * pow((pow(10, log10mu)), i_parm.kc);
            for(j=nvar; j<nvar+ncon; j++){
                k = row_strt[j]-1;
                a[k] += -i_pert->d_c;
            }
        }
        else{
            (*pert_pivot) = 1; /* Ask to change the pivot tolerance */
            i_pert->jacobian_perturbed = 1;

            fprintf(stderr, "W[K_AUG]...\t[INERTIA_STRATEGY]"
                            "d_c is already > 0\nThe Jacobian might be singular. Asking for better accuracy..\n");
        }
    }
        /*
         * else if(zi!=0){
         *         fprintf(stderr,"E[K_AUG]...\t[INERTIA_STRATEGY]"
         *                                "Failure, there is a zero eigenvalue. The jacobian is possibly singular.\n");
         *                                        exit(-1);  Would this block ever get evaluated?}*/
    else{
        i_pert->d_w_last = d_w_trial;
        i_pert->d_c_last = i_pert->d_c;
        printf("I[K_AUG]...\t[INERTIA_STRATEGY]"
               "Inertia check successful.\n");
        return 0;
    }
    (*try_n)++;
    return 1;



}

