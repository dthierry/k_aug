/* @source untitled.c
** beta 01
** June 20th, 2017
** @author: David Thierry (dmolinat@andrew.cmu) dav0@lb2016-1

********************************************************************************

@fun_name ********************************************
**
** Computes mu and adjust the primal variables value if necesary.
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

#include "mu_adjust_primal.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void mu_adjust_x(int nvar, double *x, double *lbv, real *zL, real *zU){
	/* find a multiplier-primal combination that allows to compute mu*/
	/* adjust primal */
	int i;
	double mul, muu, mu, logmu0;
	

	logmu0 = 0.0;
	for(i=0; i<nvar; i++){
		mul = 0.0;
		muu = 0.0;
		mu  = 0.0;
		if(zL[i] > 0){
			mul = (x[i] - lbv[2*i]) * (zL[i]);
		}
		if(-zU[i] > 0){
			muu = (x[i] - lbv[2*i + 1])*(zU[i]);
		}
		if(mul>0.0 && muu>0.0){
			mu = mul;
		}
		else if (mul > 0.0){
			mu = mul;
		}
		else if (muu > 0.0){
			mu = muu;
		}
		if(mu>0.0){
			if(fabs(logmu0 - log(mu)) < 1e-10){
				printf("I[KMATRIX]...\t[KMATRIX_ASL]"
					"log(mu) computed=%.g at var_i=%d\n", log(mu), i);
				break;
			}
			else{
				logmu0 = log(mu);
			}
		}
	}
	/*mu = (log(mu) > -8.6) ? exp(-8.6): mu;*/
	for(i=0; i<nvar; i++){

		if(zL[i]>0 && -zU[i] > 0){
			if((x[i] - lbv[2*i]) * (zL[i]) < mu*0.5 ||  (x[i] - lbv[2*i + 1])*(zU[i]) < mu*0.5){
				x[i] = ((x[i] - lbv[2*i]) * (zL[i])) < ((x[i] - lbv[2*i + 1])*(zU[i])) ? 
				mu/(zL[i]) + lbv[2*i]: mu/(zU[i]) +lbv[2*i+1];
			}
		}
		else if(zL[i] > 0){
			if((x[i] - lbv[2*i]) * (zL[i]) < mu*0.5){
				x[i] = mu/(zL[i]) + lbv[2*i]; 
			}
		}
		else if(-zU[i] > 0){
			if((x[i] - lbv[2*i + 1])*(zU[i]) < mu*0.5){
				x[i] = mu/(zU[i]) + lbv[2*i + 1];
			}
		}

		
	}

}

