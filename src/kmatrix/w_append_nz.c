/* @source w_append_nz.c
**
** April 5th, 2017
** @author: David Thierry (dmolinat@andrew.cmu) dav0@lb2016-1

******************************************************************************

@wnzappnd ********************************************
**
** Builds nz appended W matrix
**
** @param [r] Wrow
** @param [r] Wcol
** @param [r] Wij
** @@
*******************************************************************************/
#include "w_append_nz.h"
#include <stdio.h>
#include <assert.h>

void wnzappnd(fint *Wrow, fint *Wcol, real *Wij, fint Wnz, fint nvar,
	fint *Wr_new, fint *Wc_new, real *Wi_new, fint Wnz_new){

	int j, k, l;
	fint i;
	fint ptr0;
	FILE *fdebug;

	l = 0;
	ptr0 = 0;

	//should skip this if Wnz == 0
	for(i=1; i<=nvar; i++){
		//assert(l<Wnz_new);
		for(j=ptr0; j<Wnz; j++){
			// found col
			if(Wcol[j] == i){
				k = j;
				while(k < Wnz && Wcol[j] == Wcol[k]){
					Wr_new[l] = Wrow[k];
					Wc_new[l] = Wcol[k];
					Wi_new[l] = Wij[k];
					l++;
					k++;
					if(k == Wnz){break;}
				}
				// check main diagonal (upper triangular)
				if(Wrow[k-1] != Wcol[k-1]){
					Wr_new[l] = i;
					Wc_new[l] = i;
					Wi_new[l] = 0.0;
					l++;
				}
				if(k < Wnz){
					ptr0 = k;	
				}
				else{
					ptr0 = k-1;		
				}
				
				break;
		  }
		  // row not found appending new element
		  else if (Wcol[j] > i){
		  	Wr_new[l] = i;
				Wc_new[l] = i;
				Wi_new[l] = 0.0;
				l++;
				break;
			}
		}
	}
// Last col is less than nvar, appending missing elements
if(Wnz != 0){
	if(Wcol[Wnz-1] < nvar){
		for(i=(Wcol[Wnz-1]+1); i<=nvar; i++){
			assert(l<Wnz_new);
			Wr_new[l] = i;
			Wc_new[l] = i;
			Wi_new[l] = 0.0;
			l++;
		}
	}	
}
else{
	for(i=1; i<=nvar; i++){
			assert(l<Wnz_new);
			Wr_new[l] = i;
			Wc_new[l] = i;
			Wi_new[l] = 0.0;
			l++;
		}
}

fdebug = fopen("debug_hessian_appended.txt", "w");
for(l=0; l<Wnz_new; l++){
	fprintf(fdebug, "\t%ld\t%ld\t%.g\n", Wr_new[l], Wc_new[l], Wi_new[l]);
}
fclose(fdebug);


}
