/* @source assemble_rhsds_red_hess.c
** beta 0
** April 18th, 2017
** @author: David Thierry (dmolinat@andrew.cmu) dav0@lb2016-1

********************************************************************************

@assemble_rhsds_red_hess ********************************************
**
** Assembles right hand sides for the reduced hessian
**
** @param [r] n_rhs
** @param [r] rhs_len
** @param [r] rhsbksolv
** @param [r] dp
** @param [r] nvar
** @param [r] ncon
** @param [r] rhs_ptr
** @@
*******************************************************************************/

#include "assemble_rhs_dcdp.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

typedef struct npdp_coord {
    int row;
    int col;
} npdp_coord;

int comp_f0(const void *a, const void *b);


void assemble_rhs_dcdp(real **rhs_dcdp, fint nvar, fint ncon, int *n_p, int *n_x,
                       SufDesc *dcdp, int **hr_point, SufDesc *var_order) {
    int j, i;
    int temp;
    FILE *somefile;
    int *help_list;
    int *help_list_order;
    npdp_coord *ordered_npdp;

    /* We don't know the amount of rhs's so allocate the maximum possible namely n+m*/
    /* There is a potential error here */
    help_list = (int *) malloc(sizeof(int) * (nvar + ncon));
    help_list_order = (int *) malloc(sizeof(int) * (nvar + ncon));
    /* Then use hr_point to mark the positions (KKT) indexed by the corresponding rhs*/

    /* Making the matrix by column as usual in fotran*/

    (*n_p) = 0; /* Count number of parameters */

    for (i = 0; i < ncon; i++) {
        temp = dcdp->u.i[i]; /* Retrieve value of suffix*/
        /* Fortran-style */
        /* Find non-zero */
        if (temp != 0) {
            if ((temp - 1) > ncon) {/* error */
                printf("E[K_AUG]...\t[ASSM_RHS_DCDP]"
                       "The suffix at %d is greater than n_con e.g. %d\n", i, ncon);
                exit(-1);

            } else if ((temp - 1) < 0) {/* error again*/
                printf("E[K_AUG]...\t[ASSM_RHS_DCDP]"
                       "The suffix at %d is negative (%d)\n", i, temp);
                exit(-1);

            }
            help_list_order[(*n_p)] = temp; /* Save the position of nz */
            help_list[(*n_p)] = i + nvar; /* Save the position of nz */
            (*n_p)++;
        }
    }

    printf("I[K_AUG]...\t[ASSM_RHS_DCDP]"
           "According to the suffixes declared len p is %d \n", *(n_p));
    if ((*n_p) <= (0)) {
        printf("E[K_AUG]...\t[ASSM_RHS_DCDP]"
               "No valid number of n_p declared\n");
        exit(-1);
    }

    ordered_npdp = (npdp_coord *) malloc(sizeof(npdp_coord) * (*n_p));
    for (i = 0; i < (*n_p); i++) {
        (ordered_npdp + i)->row = help_list[i];
        (ordered_npdp + i)->col = help_list_order[i];
    }

    /*for(i=0; i<(*n_p); i++){printf("%d\t%d\t%d\n", i, (ordered_npdp+i)->row, (ordered_npdp+i)->col);}*/
    qsort(ordered_npdp, (*n_p), sizeof(ordered_npdp), comp_f0);
    /*printf("Let us do it again..\n");*/
    /*for(i=0; i<(*n_p); i++){printf("%d\t%d\t%d\n", i, (ordered_npdp+i)->row, (ordered_npdp+i)->col);}*/


    (*rhs_dcdp) = (real *) calloc((nvar + ncon) * (*n_p), sizeof(real));
    (*hr_point) = (int *) malloc(sizeof(int) * (*n_p));

    memset((*rhs_dcdp), 0, (nvar + ncon) * (*n_p) * sizeof(real));
    /* Put a 1.0 in the appropiate position.*/
    for (i = 0; i < (*n_p); i++) {
        /*(*hr_point)[i] = help_list[i];*/
        (*hr_point)[i] = (ordered_npdp + i)->row;
        /*printf("i %d help %d\n", i, help_list[i]);*/
        (*rhs_dcdp)[(nvar + ncon) * i + (*hr_point)[i]] = -1.0;
    }

    free(help_list);
    free(help_list_order);


    somefile = fopen("rhs_dcdp", "w");


    for (i = 0; i < (*n_p); i++) {
        fprintf(somefile, "\t%d\t", (*hr_point)[i]);
    }
    fprintf(somefile, "\n\n");

    for (i = 0; i < (nvar + ncon); i++) {
        for (j = 0; j < (*n_p); j++) {
            fprintf(somefile, "\t%f\t", *((*rhs_dcdp) + j * (nvar + ncon) + i));
        }
        fprintf(somefile, "\n");
    }
    fclose(somefile);



    /* order by row, if the var_order suffix was declared*/
    if (var_order->u.i) {
        j = 0;
        for (i = 0; i < nvar; i++) { if (var_order->u.i[i] > 0) { j++; }}
        (*n_x) = j; /* number of vars of interest */
        ordered_npdp = (npdp_coord *) realloc(ordered_npdp, sizeof(npdp_coord) * (*n_x));
        j = 0;
        for (i = 0; i < nvar; i++) {
            if (var_order->u.i[i] > 0) {
                (ordered_npdp + j)->col = var_order->u.i[i];
                (ordered_npdp + j)->row = i;
                j++;
                /*if(j>=(*n_x)){break;}*/
            }
        }
    }

    qsort(ordered_npdp, (*n_x), sizeof(ordered_npdp), comp_f0);
    (*hr_point) = (int *) realloc((*hr_point), sizeof(int) * (*n_x));
    memset((*hr_point), 0, sizeof(int) * (*n_x));
    for (i = 0; i < (*n_x); i++) { (*hr_point)[i] = (ordered_npdp + i)->row; }
    /* the order will go to hr_point */

    free(ordered_npdp);
}


int comp_f0(const void *a, const void *b) {
    npdp_coord atemp = *(const npdp_coord *) a;
    npdp_coord btemp = *(const npdp_coord *) b;
    int aint = (int) (atemp.col);
    int bint = (int) (btemp.col);
    return (aint - bint);
}
