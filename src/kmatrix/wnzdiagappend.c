#include "wnzdiagappend.h"
#include <stdio.h>
#include <assert.h>

void wnzdappend(fint *Wrow, fint *Wcol, real *Wij, fint Wnz, fint Wnrow);


void wnzdappend(fint *Wrow, fint *Wcol, real *Wij, fint Wnz, fint Wnrow){
	// the number of missing non-zeroes
	int nMz;
	fint row=1;
	fint nZmd=0;
	fint newNZ=0;
	fint *Wrowa, *Wcola;
	real *Wija;
	fint i;
	int j, k, l;
	unsigned char flag;
	fint ptr0;
	ptr0 = 0;
	int miss_nz, miss_row;
	miss_row = 0;
	miss_nz = 0;
	FILE *fdebug;

	// count the number of missing nz
	// first loop on vector elements (starts 1)
	for(i=1; i < Wnrow+1; i++){
		for(j=ptr0; j<Wnz ;j++){
			// found row
			if(Wrow[j] == i){
				k = j;
				while(k < Wnz-1 && Wrow[k] == Wrow[j]){
					// found nz in main diag
					if(Wcol[k] == i){
						break;
					}
					k++;
				}
				ptr0 = k;
				if(Wcol[k] == i){
					break;
				}
				// not found nz in main diag
				else{
					miss_nz++;
					break;}
			}
			// not found row
			else if (Wrow[j] > i){
				miss_row++;
				break;
			}
		}
	}

	printf("Number of elements missing in the main diag %d\n", miss_nz);
	printf("Number of rows missing %d\n", miss_row);
	newNZ = miss_row + miss_nz + Wnz;

	// allocate new array
	Wrowa = (fint *)malloc(sizeof(fint)*(newNZ));
	Wcola = (fint *)malloc(sizeof(fint)*(newNZ));
	Wija = (real *)malloc(sizeof(real)*(newNZ));
	l = 0;
	ptr0 = 0;
	for(i=1; i<=Wnrow; i++){
		assert(l<newNZ);
		for(j=ptr0; j<Wnz; j++){
			// found row
			if(Wrow[j] == i){
				// check main diagonal (upper triangular)
				if(Wrow[j] == Wcol[j]){
					Wrowa[l] = Wrow[j];
					Wcola[l] = Wcol[j];
					Wija[l] = Wij[j];
					l++;
				}
				else{
					Wrowa[l] = i;
					Wcola[l] = i;
					Wija[l] = 0.0;
					l++;
				}
				k = j + 1;
				// update rest of the matrix
				while(Wrow[k] == Wrow[j] && k < Wnz-1){
					Wrowa[l] = Wrow[k];
					Wcola[l] = Wcol[k];
					Wija[l] = Wij[k];
					k++;
					l++;	
				}
				ptr0 = k;
				break;
		  }
		  // row not found appending new element
		  else if (Wrow[j] > i){
		  	Wrowa[l] = i;
				Wcola[l] = i;
				Wija[l] = 0.0;
				l++;
				break;
			}
		}
	}

	fdebug = fopen('debug_hessian_appended.txt', 'w');
	for(l=0; l<newNZ; l++){
		fprintf(fdebug, "\t%ld\t%ld\t%f\n", Wrowa[l], Wcola[l], Wija[l]);
	}
	fclose(fdebug);

	free(Wrowa);
	free(Wcola);
	free(Wija);
	/*
	flag = 0;
	for(int i=0; i<Wnz; i++){
		j = 0;
		k = 0;
		if (Wrow[i] == row){
			if(flag){
				// off-diagonal
				Wrowa[j] = Wrow[i];
				Wcola[j] = Wcol[i];
				Wija[j] = Wij[i];
				j++;
				continue;
			}
			if(Wrow[i] == Wcol[i]){
				// element on main diagonal
				Wrowa[j] = Wrow[i];
				Wcola[j] = Wcol[i];
				Wija[j] = Wij[i];
				flag = 1;
			}
			else{
				// no element on main diagonal
				Wrowa[j] = row;
				Wcola[j] = row;
				Wija[j] = 0.0;
				flag = 1;
				j++;
				Wrowa[j] = Wrow[i];
				Wcola[j] = Wcol[i];
				Wija[j] = Wij[i];
			}
			j++;

			k++;
		}
		else{
			flag = 0;
			if(Wrow[i]-row > 2){
				printf("Warning row skip of %ld\n", Wrow[i]-row);
				printf("Where %ld\n", Wrow[i]);
			}
			row++;
			if(Wrow[i] == Wcol[i]){
				// do smth
				Wrowa[j] = Wrow[i];
				Wcola[j] = Wcol[i];
					Wija[j] = Wij[i];
					flag = 1;
				}
			else{
				Wrowa[j] = row;
				Wcola[j] = row;
				Wija[j] = 0.0;
				flag = 1;
				j++;
				Wrowa[j] = Wrow[i];
				Wcola[j] = Wcol[i];
				Wija[j] = Wij[i];
			}
		j++;
		}
	}

	*/
}