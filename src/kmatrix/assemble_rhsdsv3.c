/* @source assemble_rhsdsv3.c
** beta 0
** April 18th, 2017
** @author: David Thierry (dmolinat@andrew.cmu) dav0@lb2016-1

********************************************************************************

@assemble_rhsds ********************************************
**
** Assembles right hand sides from the suffixes in asl. Also assembles delta p
** The right hand sides have an identity matrix as the derivative of c(x;p) wrt p
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

#include "assemble_rhsdsv3.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

void assemble_rhsds(int n_rhs, fint rhs_len, 
 real *rhsbksolv, real *dp, fint nvar, fint ncon, SufDesc **rhs_ptr){
 	int j, k, i;
 	real *temp;
 	FILE *somefile;
 	/* by [n][rhs] */
 	/* Firts n elements are zero */
 	/*for(j=0; j < n_rhs; j++){ */
 	/*	for(k=0; k < nvar; k++){ */
 	/*		temp = *(rhsbksolv + k); */
 	/*		temp[j] = 0.0; */
 	/*	} */
 	/*} */
 	
 	for(j=0; j < n_rhs; j++){
 		i = 0;
    /*incidentally storing the matrix by column, as fortran demands*/
 		for(k=nvar; k < (nvar+ncon); k++){
 			temp = (rhsbksolv + k + rhs_len*j);
 			/* *temp = 0.0; ((*(rhs_ptr+j))->u.r[k-nvar]); */ 
                  if ( ((*(rhs_ptr+j))->u.r[k-nvar]) != 0.0 ){
                        *temp = -1.0; /* this must be negative*/
                  			*(dp + j) = ((*(rhs_ptr+j))->u.r[k-nvar]);
                        printf("rhs_%d, Gc_P at %d\n", j, k-nvar);
                        i++;
                  }
    if (i > 1){
    	printf("W[KMATRIX]...\t[ASSEMBLE_RHSDS]"
        "Current rhs_%d has more than one element!! \n" , j);
    }
 			/**temp = ((*(rhs_ptr+j))->u.r[k-nvar]);*/
 			/*printf("%f\n", *temp);*/
 			/*assert( (*(rhs_ptr+j))->u.r[k-nvar] == 0.0 );*/
 		}
 		
 	}
 	

 	somefile = fopen("rhs_sens", "w");

 	for(j=0; j < (nvar+ncon); j++){
 		for(k=0; k < n_rhs; k++){
 			fprintf(somefile, "\t%f\t", *(rhsbksolv + j + rhs_len*k) );
 		}
 		fprintf(somefile, "\n");
 	}
 	fclose(somefile);


}
