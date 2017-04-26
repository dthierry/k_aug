#include "assemble_rhsds.h"
#include <stdio.h>
#include <stdlib.h>

void assemble_rhsds(int n_rhs, fint rhs_len, 
 real **rhsbksolv, fint nvar, fint ncon, SufDesc *rhs_ptr){
 	int j, k, l;
 	real *temp;
 	l = 0;

 	for(j=0; j < nvar; j++){
 		for(k=0; k < n_rhs; k++){
 			temp = *(rhsbksolv + k);
 			temp[j] = 0.0;
 		}
 	}
 	
 	for(j=nvar; j < (nvar+ncon); j++){
 		for(k=0; k < n_rhs; k++){
 			temp = *(rhsbksolv + k);
 			temp[j] = (real)((rhs_ptr+k)->u.r[l]);
 		}
 		l++;
 	}
 	
 	for(j=nvar; j < (nvar+ncon); j++){

 		printf("%f\t\n",rhsbksolv[0][j] );
 	}


 	//for(j=0; j < ncon; j++){
  //temp = rhs_ptr->u.i[j];
 	//printf("rhs 0 %d\n", temp);
 	//}
	
}