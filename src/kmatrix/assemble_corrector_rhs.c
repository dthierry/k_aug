#include "assemble_corrector_rhs.h"
#include "get_grad_f.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

// must traverse jac matrix by row
/* @source k_assemble.c
**
** April 25th, 2017
** @author: David Thierry (dmolinat@andrew.cmu) dav0@lb2016-1

******************************************************************************

@assemble_corrector_rhs ********************************************
**
** Assembles [df + lamda*dC; cA] (dense)
**
** @param [r] Wrow
** @param [r] Wcol
** @param [r] Wij
** @@
*******************************************************************************/

void assemble_corrector_rhs(ASL *asl, real *x, real *lambda,
  fint nvar, fint ncon,
	fint *Arow, fint *Acol, real *Aij, fint Anz,
  real **Crhs, fint *nerror, int *c_flag){
	
	fint i;
	int j, k, l, ptr;
	// column starts
	int *cs_a;
	// gradient part
	real *rhs;
	// constraint part
	real *rhs_con;
	// full rhs
	real *rhs_full;

	// A row and col transposes
	fint *Ar_t, *Ac_t;

	FILE *somefile;


	cs_a = (int *)malloc(sizeof(int) * ncon);
	assert(cs_a != NULL);

	// We want A but ASL gives us J, so we have to transpose;
	Ar_t = Acol;
	Ac_t = Arow;


	// Initialize the row starts with -1 just in case;

	for(j=0; j<ncon; j++){cs_a[j] = -1;}

	// traverse A list to find col starts
	ptr = 0;
	for(i=1; i<=ncon; i++){
		for(k=ptr; k<Anz; k++){
			// found col
			if(Ac_t[k]==i){
				cs_a[i-1] = k;
				ptr = k;
				break;
			}
		}
	}

	// get gradient of f
	rhs = (real *)malloc(sizeof(real) * (nvar));
	


	objgrd(0, x, rhs, nerror);
	
	

	rhs_con = (real *)malloc(sizeof(real) * (ncon));

	conval(x, rhs_con, nerror);


	// multiply by lambda(s)
	for(i=1; i<=nvar; i++){
		// Traverse the A matrix by col;
		for(j=1; j<=ncon; j++){
			k = cs_a[j-1];
			// if there is no column
			if(k == -1){continue;}
			l = k;
			// look for element on the column
			while(Ac_t[k]==Ac_t[l] && l < Anz){
				// found element in column
				if(Ar_t[l]==i){
					// move the current pointer so we dont start from the beggining
					// every time
					cs_a[j-1] = l;
					// add dC  
					rhs[i-1] += Aij[l] * lambda[j-1];
					break;}
				l++;
				// must put this here otw there will be a read beyond
				// array
				if(l == Anz){break;}
			}
		}
	}


//	for(i = 0; i<(nvar); i++){
		//printf("irhs %d, val= %f \n",i , rhs[i]);
	//}

	//for(i = 0; i<(ncon); i++){
		//printf("irhs %d, val= %f \n",i , rhs_con[i]);
	//}
	rhs_full = (real *)calloc((ncon + nvar), sizeof(real));	

	for(i = 0; i<(nvar); i++){
		rhs_full[i] = rhs[i];
	}
	j = 0;
	for(i = nvar; i<(nvar + ncon); i++){
		// rhs_full[i] = rhs_con[j];
		rhs_full[i] = -rhs_con[j];
		if(c_flag[j] == -3){
			// rhs_full[i] += -LUrhs[2*j]; original
			rhs_full[i] += LUrhs[2*j];
		}
		else if(c_flag[j] == -1){
			// rhs_full[i] += -LUrhs[2*j+1];
			rhs_full[i] += LUrhs[2*j+1];
		}
		else if(c_flag[j] == -2){
			// rhs_full[i] += -LUrhs[2*j];
			rhs_full[i] += LUrhs[2*j];
		}
		
		j++;
	}
	somefile = fopen("rhs_corrector.txt", "w");

	for(i=0; i<(nvar + ncon); i++){
		fprintf(somefile, "%d\t%f\n", i, rhs_full[i]);
	}
	fclose(somefile);

	free(cs_a);
	free(rhs);
	free(rhs_con);
	*Crhs = rhs_full;
}
