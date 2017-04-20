#include "find_inequalities.h"
#include <stdio.h>
#include <stdlib.h>


void find_ineq_con(fint ncon_,real *LBC){
	int i;

	for(i=0; i<ncon_; i++){
		printf("%d LBC %f \t UBC %f\n",i, LBC[2*i], LBC[2*i+1]);
	}
}



void find_bounds(fint nvar_,real *lbv){
	int i;

	for(i=0; i<nvar_; i++){
		printf("%d LBC %f \t UBC %f\n",i, lbv[2*i], lbv[2*i+1]);
	}
}
