#include "wcreord.h"
#include <stdio.h>
#include <assert.h>
void w_col_sort(fint *Wrow, fint *Wcol, real *Wij, fint Wnz, fint Wnrow){
	fint i;
	int j, k, l, m, ptr;
	int cs[Wnrow];
	ptr = 0;
	fint Wrt[Wnz], Wct[Wnz];
	FILE *somefile, *cprint;
	real Wijt[Wnz];

	for(j=0; j<Wnrow; j++){cs[j] = -1;}

	// traverse list to find col starts
	for(i=1; i<=Wnrow; i++){
		if(i == 258){
			printf("Stop\n");
		}
		for(k=ptr; k<Wnz; k++){
			// found col
			if(Wcol[k]==i){
				cs[i-1] = k;
				ptr = k;
				break;
			}
		}
	}

	cprint = fopen("c_.txt", "w");
	for(i=1; i<=Wnrow; i++){
		fprintf(cprint, "%d\n", cs[i-1]);
	}
	fclose(cprint);

	// reorder list
	m = 0;
	for(i=1; i<=Wnrow; i++){
		if(i==33){
			printf("Something\n");
		}
		for(j=1; j<=Wnrow; j++){
			k = cs[j-1];
			// if there is no column
			if(k == -1){continue;}
			l = k;
			// look for element on the column
			while(Wcol[k]==Wcol[l] && l< Wnz && i<=Wcol[l]){
				// found element in column
				if(Wrow[l]==i){
					cs[j-1] = l;
					Wrt[m] = Wrow[l];
					Wct[m] = Wcol[l];
					Wijt[m] = Wij[l];
					m++;
					break;}
				l++;
				// must put this here otw there will be a read beyond
				// array
				if(l == Wnz){break;}
			}
		}
	}
	somefile = fopen("reordered_hessian.txt", "w");
	
	for(j=0; j<Wnz; j++){fprintf(somefile, "\t%ld\t%ld\t%.g\n", Wrt[j], Wct[j], Wijt[j]);}
	fclose(somefile);
	
}