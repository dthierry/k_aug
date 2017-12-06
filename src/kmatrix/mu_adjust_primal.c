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
	/*If zL and zU are not declared, they will have a 0.0 value, presumably
	  therefore mul and muu will have a 0.0 value as well.
	*/

	logmu0 = -11.0;
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
			if(fabs(logmu0 - log(mu)) < 1){
				printf("I[KMATRIX]...\t[ADJUST_MU]"
					"log(mu) close to the target\t%.g at var_i=%d\n", log(mu), i);
                        logmu0 = log(mu);
				break;
			}
			else{
				logmu0 = log(mu);
				printf("I[KMATRIX]...\t[ADJUST_MU]"
					"log(mu) computed=%.g at var_i=%d\n", log(mu), i);
			}
            }
            else if ((fabs(x[i] - lbv[2*i]) < 1e-08) && (fabs(zL[i]) < 1e-08)){
                  printf("I[KMATRIX]...\t[ADJUST_MU]"
		      	 "\tWarning strict complementarity (lb) is not within tol for var:%d\n", i);

            }
            else if ((fabs(x[i] - lbv[2*i + 1]) < 1e-08) && (fabs(zU[i]) < 1e-08)){
                  printf("I[KMATRIX]...\t[ADJUST_MU]"
		      	 "\tWarning strict complementarity (ub) is not within tol for var:%d\n", i);

            }


	}
      if(logmu0 == -11.0){
            printf("I[KMATRIX]...\t[ADJUST_MU]"
			 "\tWarning no relevant info from the problem can predict logmu\n");
      }
      else if(logmu0 > -11.0){
            printf("I[KMATRIX]...\t[ADJUST_MU]"
			 "\tWarning logmu is over the target; is this optimal?\n");
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

