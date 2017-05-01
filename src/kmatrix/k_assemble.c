/* @source k_assemble.c
**
** April 5th, 2017
** @author: David Thierry (dmolinat@andrew.cmu) dav0@lb2016-1

******************************************************************************

@k_assemble ********************************************
**
** Assembles [W A 0]
**
** @param [r] Wrow
** @param [r] Wcol
** @param [r] Wij
** @@
*******************************************************************************/

#include "k_assemble.h"
#include <assert.h>

void k_assemble(fint *Wrow, fint *Wcol, real *Wij, fint Wnz, fint nvar, fint ncon,
 fint *Arow, fint *Acol, real *Aij, fint Anz,
  fint *Krow, fint *Kcol, real *Kij, fint Knz){
	
	fint i;
	int j, k, l, m, ptr;
	// column starts
	int cs_w[nvar];
	int cs_a[ncon];
	
	FILE *somefile, *cprint;

	// A row and col transposes
	fint *Ar_t, *Ac_t;

	// We want A but ASL gives us J, so we have to transpose;
	Ar_t = Acol;
	Ac_t = Arow;

	// Initialize the row starts with -1 just in case;
	for(j=0; j<nvar; j++){cs_w[j] = -1;}
	for(j=0; j<ncon; j++){cs_a[j] = -1;}

	// traverse w list to find col starts
	ptr = 0;
	for(i=1; i<=nvar; i++){
		for(k=ptr; k<Wnz; k++){
			// found col
			if(Wcol[k]==i){
				cs_w[i-1] = k;
				ptr = k;
				break;
			}
		}
	}

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


	cprint = fopen("c_.txt", "w");
	fprintf(cprint, "Hessian col starts\n");
	for(i=1; i<=nvar; i++){
		fprintf(cprint, "%d\n", cs_w[i-1]);
	}
	fprintf(cprint, "A col starts\n");
	for(i=1; i<=ncon; i++){
		fprintf(cprint, "%d\n", cs_a[i-1]);
	}
	fclose(cprint);

	// reorder list(s)
	m = 0;
	
	for(i=1; i<=nvar; i++){
		for(j=1; j<=nvar; j++){
			k = cs_w[j-1];
			// if there is no column
			if(k == -1){continue;}
			l = k;
			// look for element on the column
			while(Wcol[k]==Wcol[l] && l< Wnz && i<=Wcol[l]){
				// found element in column
				if(Wrow[l]==i){
					cs_w[j-1] = l;
					Krow[m] = Wrow[l];
					Kcol[m] = Wcol[l];
					Kij[m] = Wij[l];
					m++;
					break;}
				l++;
				// must put this here otw there will be a read beyond
				// array
				if(l == Wnz){break;}
			}
		}

		// having traversed the whole W matrix(row), proceed to traverse the A matrix by row;
		for(j=1; j<=ncon; j++){
			k = cs_a[j-1];
			// if there is no column
			if(k == -1){continue;}
			l = k;
			// look for element on the column
			while(Ac_t[k]==Ac_t[l] && l < Anz){
				// found element in column
				if(Ar_t[l]==i){
					cs_a[j-1] = l;
					Krow[m]  = Ar_t[l];
					Kcol[m]  = Ac_t[l] + ncon;
					Kij[m] = Aij[l];
					m++;
					break;}
				l++;
				// must put this here otw there will be a read beyond
				// array
				if(l == Anz){break;}
			}
		}
		//assert(m<(Knz-ncon));
	}
	// 0 part, main diagonal
	for(i=1; i<=ncon; i++){
		Krow[m] = nvar + i;
		Kcol[m] = nvar + i;
		Kij[m] = 0.0;		
		m++;
		assert(m<=Knz);
	}
	somefile = fopen("KKT.txt", "w");
	for(j=0; j<Knz; j++){fprintf(somefile, "\t%ld\t%ld\t%.g\n",
	Krow[j], Kcol[j], Kij[j]);}
	fclose(somefile);
	
}