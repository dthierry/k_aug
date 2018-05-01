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
#define HUGE_NUMBER 1e300


void compute_sigma(ASL *asl, fint nvar, real *x, SufDesc *suf_zL, SufDesc *suf_zU, real *z_L, real *z_U, real *sigma, real not_zero){
	int i;
	real sl, su;
	FILE *s_file;
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
		else{
			sl = z_L[i]/(x[i] - LUv[2*i]);
		}
		if(LUv[2*i + 1] > HUGE_NUMBER){
			su = 0.0;
		}
		else{
			/* now with correct expression */
			su = z_U[i]/(x[i] - LUv[2*i+1]);
		}
		sigma[i] = sl + su;
	}
	s_file = fopen("sigma_out", "w");
	fprintf(s_file, "%s\t%s\t%s\t%s\t%s\t%s\t\t%s\n", "i", "z_L[i]", "z_U[i]", "LUv[2*i]", "x[i]", "LUv[2*i+1]", "sigma[i]");
	for(i=0; i<nvar; i++){
		fprintf(s_file, "%d\t%.g\t%.g\t%.g\t%.g\t%.g\t\t%e %c\n", i, z_L[i], z_U[i], LUv[2*i], x[i], LUv[2*i+1], sigma[i], sigma[i]>1e-02 ? '*':' ');
	}
	fclose(s_file);

}

