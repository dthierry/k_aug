/*THIS HAS TO DO THE FOLLOWING: 
	REORDER GRADIENTS ** done
	GET NUMBER OF NZ BY ROW ** done
	NEED TO PASS NVAR & NCON ** done
	*/
#include "config_kaug.h"
#include "get_jac_asl_aug.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>


int compf_v2(const void *t1, const void *t2)
{
    temp_v2 firstt  = *(const temp_v2 *)t1;
    temp_v2 secondt = *(const temp_v2 *)t2;
    int fst = (int) (firstt.c);
    int scd = (int) (secondt.c);
    return (fst - scd);
}



int get_jac_asl_aug(ASL *asl, real *x, fint *Acol, fint *Arow, real *Aij,
                     int nvar, int ncon, fint nzc_, fint *nerror, int **nz_row_a){
    int i, j, k, l;
    cgrad *cg, **cgx;
    FILE *f_jac, *somefile;
    real *Jcont; /* container for Jac */
    temp_v2 *ttemp_v2;

    (*nz_row_a) =     (int *)calloc(sizeof(int), nvar);
    assert((*nz_row_a) != NULL);
    memset((*nz_row_a), 0, sizeof(int) * nvar);

    if(nzc_ <= 0){
        printf("E[K_AUG]...\t[GET_JAC_ASL]"
               "The Jacobian has no structural non-zeroes; exiting\n");
        /*exit(-1);*/
        return 1;
    }

    ttemp_v2 = (temp_v2 *)calloc(sizeof(temp_v2), nvar);
    Jcont =       (real *)malloc(sizeof(real) * nzc_);

    assert(Jcont != NULL);

    j = 0;
    /* Jacobian */
    /* check if there are non-zeroes in the Jacobian */


    /* evaluates the Jacobian matrix */
    jacval(x, Jcont, nerror);

    if(*nerror != 0){
        printf("E[K_AUG]...\t[GET_JAC_ASL]nerror points to %ld\n", *nerror);
        exit(-1);
    }
#ifndef PRINT_VERBOSE
    f_jac = fopen("jacobi_debug.in", "w");
    /* for analysis with mc58 */
    fprintf(f_jac, "%d\t%d\t%d\n", ncon, nvar, nzc_);
#endif

    cgx = Cgrad; /* point to the beggining of the list */
    /* Cgrad is an array of pointers
        it contains m pointers (1 for each constraint)
        each pointer defines a linked-list of gradients */
    /* We need to sort the gradient for each constraint; that is why we use
    qsort*/
    for(i = 0; i < ncon; i++) {
        /* moves by constraint */
        if(cgx == NULL){
            fprintf(stderr, "E[K_AUG]...\t[GET_JAC_ASL]"
                            "NULL GRADIENT DETECTED CONSTRAINT %d\n", i);
            exit(-1);
        }
        k = 0;
        for(cg = cgx[i]; cg; cg=cg->next){
            /* moves by nz in the constraint */
            /* actual jacobian (instead of gradient) */
#ifndef PRINT_VERBOSE
            fprintf(f_jac, "%d\t%d\t%.g\n", i+1, cg->varno+1, Jcont[cg->goff]);
#endif
            Arow[j] = cg->varno+1;
            Acol[j] = i+1;
            Aij[j] = Jcont[cg->goff];



            /* structure for qsort */
            ttemp_v2[k].c = *(Arow+j);
            ttemp_v2[k].a = *(Aij +j);
            k++;
            /* */

            (*nz_row_a)[cg->varno]++; /* Number of nz in row */
            j++;
        }

        /* ADD SMTH HERE IF U WANT INEQUALITIES */

        if(k > 1){
            /*
            for(l=0; l<k; l++){
                printf("ttemp_v2 i,l %d, %d, %d, %f\n",i+1, l, (ttemp_v2[l].c), (ttemp_v2[l].a));
            }
            */
            qsort(ttemp_v2, k, sizeof(temp_v2), compf_v2); /* Do the sorting.
	  	This actually sorts the array of temp_v2s*/
            for(l=0; l<k; l++){
                /* Re-write the current column */
                Arow[(j-1) - (k-1) + l] = (ttemp_v2[l].c);
                Aij [(j-1) - (k-1) + l] = (ttemp_v2[l].a);
                /*
                printf("ttemp_v2 i,l %d, %d, %d, %f\n",i+1, l, (ttemp_v2[l].c), (ttemp_v2[l].a));
                */
            }
            memset(ttemp_v2, 0, sizeof(temp_v2)* k);
        }
    }
#ifndef PRINT_VERBOSE
    fclose(f_jac);
#endif
#ifndef PRINT_VERBOSE
    somefile = fopen("anz_p_row.in", "w");
    for(i=0; i<nvar; i++){fprintf(somefile, "%d\n", (*nz_row_a)[i]);}
    fclose(somefile);
#endif
    free(Jcont);
    free(ttemp_v2);
#ifndef PRINT_VERBOSE
    f_jac = fopen("grad_debug_sorted.in", "w");
    for(i=0; i<j; i++){
        fprintf(f_jac, "%d\t%d\t%.g\n", Arow[i] , Acol[i], Aij[i]);
    }
    fclose(f_jac);
#endif
    return 0;

}
