#include "wapnz.h"
#include <stdio.h>
#include <assert.h>
#include "wcreord.h"


void wnzappnd(fint *Wrow, fint *Wcol, real *Wij, fint Wnz, fint Wnrow){
	// the number of missing non-zeroes
	//int nMz;
	//fint row=1;
	//fint nZmd=0;
	fint newNZ=0;
	fint *Wrowa, *Wcola;
	real *Wija;
	fint i;
	int j, k, l;
	//unsigned char flag;
	fint ptr0;
	ptr0 = 0;
	int miss_nz, miss_row;
	miss_row = 0;
	miss_nz = 0;
	FILE *fdebug;

	// count the number of missing nz
	// first loop on vector elements (starts 1)
	for(i=1; i <= Wnrow; i++){
		for(j=ptr0; j<Wnz ;j++){
			// found col
			if(Wcol[j] == i){
				k = j;
				while(k < Wnz-1 && Wcol[k] == Wcol[j]){
					// found nz in main diag
					if(Wrow[k] == i){
						break;
					}
					k++;
				}
				ptr0 = k;
				if(Wrow[k] == i){
					break;
				}
				// not found nz in main diag
				else{
					miss_nz++;
					break;
				}
			}
			// not found col
			else if (Wcol[j] > i){
				miss_row++;
				break;
			}
			//else{
				//printf("Something else happened\n");
			//}
		}
	}

	printf("Number of elements missing in the main diag %d\n", miss_nz);
	printf("Number of rows missing %d\n", miss_row);
	newNZ = miss_row + miss_nz + Wnz;

	if(Wcol[Wnz-1] < Wnrow){
		printf("Last element provided is less the square matrix nCols\n");
		printf("Offset of %ld\n", Wnrow - Wcol[Wnz-1]);
	}

	// allocate new array
	Wrowa = (fint *)malloc(sizeof(fint)*(newNZ));
	Wcola = (fint *)malloc(sizeof(fint)*(newNZ));
	Wija = (real *)malloc(sizeof(real)*(newNZ));
	l = 0;
	ptr0 = 0;
	for(i=1; i<=Wnrow; i++){
		//assert(l<newNZ);
		if(i == 1517){
					printf("here\n");
			}
		for(j=ptr0; j<Wnz; j++){
			// found col
			if(Wcol[j] == i){
				k = j;
				while(k < Wnz && Wcol[j] == Wcol[k]){
					Wrowa[l] = Wrow[k];
					Wcola[l] = Wcol[k];
					Wija[l] = Wij[k];
					l++;
					k++;
					if(k == Wnz){break;}
				}
				// check main diagonal (upper triangular)
				if(Wrow[k-1] != Wcol[k-1]){
					Wrowa[l] = i;
					Wcola[l] = i;
					Wija[l] = 0.0;
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
		  	Wrowa[l] = i;
				Wcola[l] = i;
				Wija[l] = 0.0;
				l++;
				break;
			}
		}
	}

fdebug = fopen("debug_hessian_appended.txt", "w");
for(l=0; l<newNZ; l++){
	fprintf(fdebug, "\t%ld\t%ld\t%.g\n", Wrowa[l], Wcola[l], Wija[l]);
}
fclose(fdebug);

w_col_sort(Wrowa, Wcola, Wija, newNZ, Wnrow);


free(Wrowa);
free(Wcola);
free(Wija);

}