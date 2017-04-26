#include "find_inequalities.h"
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <assert.h>

void find_ineq_con(fint ncon_,real *LBC, int *c_flag){
	int i;
	int geq_con, leq_con, eq_con;
	geq_con=0;
	leq_con=0;
	eq_con =0;
	for(i=0; i<ncon_; i++){
		//printf("%d LBC %f \t UBC %f\n",i, LBC[2*i], LBC[2*i+1]);
		//assert(LBC[2*i] <= -1e300);
		if(LBC[2*i] <= -1e300){
			leq_con++;
			c_flag[i] = 1;
			if(LBC[2*i+1]!=0.0){c_flag[i] = -1;}
		}
		else if(LBC[2*i+1] >= 1e300){
			geq_con++;
			c_flag[i] = 2;
			if(LBC[2*i]!=0.0){c_flag[i] = -2;}
		}
		else if (LBC[2*i] <= -1e300 && LBC[2*i+1] >= 1e300){
			printf("I[KMATRIX]...\t[FIND_INEQUALITIES]"
				"-inf <= c(x) <= +inf detected\b");
			c_flag[i] = 299;
		}
		else{
			eq_con++;
			c_flag[i] = 3;
			if(LBC[2*i]!=0.0){c_flag[i] = -3;}
		}
	}
	printf("I[KMATRIX]...\t[FIND_INEQUALITIES]"
				"summary: eq: %d, leq: %d, geq: %d \n", eq_con, leq_con, geq_con);
}



void find_bounds(fint nvar_,real *lbv){
	int i;

	for(i=0; i<nvar_; i++){
		printf("%d LBC %f \t UBC %f\n",i, lbv[2*i], lbv[2*i+1]);
	}
}
