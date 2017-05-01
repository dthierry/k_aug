#include "assemble_rhsds.h"
#include <stdio.h>
#include <stdlib.h>

void assemble_rhsds(int n_rhs, fint rhs_len, 
 real **rhsbksolv, fint nvar, fint ncon, SufDesc **rhs_ptr){
 	int j, k, l;
 	real *temp;
 	FILE *somefile;
 	// by [n][rhs]

 	//for(j=0; j < n_rhs; j++){
 	//	for(k=0; k < nvar; k++){
 	//		temp = *(rhsbksolv + k);
 	//		temp[j] = 0.0;
 	//	}
 	//}
 	
 	for(j=0; j < n_rhs; j++){
 		for(k=nvar; k < (nvar+ncon); k++){
 			temp = *(rhsbksolv + k);
 			temp[j] = ((*(rhs_ptr+j))->u.r[k-nvar]);
 		}
 		
 	}
 	

 	somefile = fopen("rhs_sens", "w");

 	for(j=0; j < (nvar+ncon); j++){
 		for(k=0; k < n_rhs; k++){
 			fprintf(somefile, "\t%f\t", *(*(rhsbksolv + j)+k) );
 		}
 		fprintf(somefile, "\n");
 	}
 	fclose(somefile);

 	//for(j=0; j < ncon; j++){
  //temp = rhs_ptr->u.i[j];
 	//printf("rhs 0 %d\n", temp);
 	//}
	
}
