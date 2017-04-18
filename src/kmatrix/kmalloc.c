/* @source KKT_matrix
**
** April 5th, 2017
** @author: David Thierry (dmolinat@andrew.cmu) dav0@lb2016-1

******************************************************************************

@k_malloc tcode_readdata ********************************************
**
** Read Calculate space for KKT matrix
**
** @param [r] Wrow
** @param [r] Wcol
** @return [fint] Space required
** @@
*******************************************************************************/
#include "kmalloc.h"
// Calculates the amount of space for the whole KKT matrix
fint k_malloc(fint *Wrow, fint *Wcol, real *Wij, fint Wnz, fint nvar, fint ncon,
 fint Anz){

	fint newNZ;
	fint i;
	int j, k;

	//unsigned char flag;
	fint ptr0;
	ptr0 = 0;
	int miss_nz, miss_row;
	
	miss_row = 0;
	miss_nz = 0;
	fint Knz;

	// count the number of missing nz
	// first loop on vector elements (starts 1)
	for(i=1; i <= nvar; i++){
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

	printf("Number of elements missing in the main diag W%d\n", miss_nz);
	printf("Number of rows missing W %d\n", miss_row);
	newNZ = miss_row + miss_nz + Wnz;

	if(Wcol[Wnz-1] < nvar){
		printf("Last element provided is less the square matrix nCols W\n");
		printf("Offset of %ld\n", nvar - Wcol[Wnz-1]);
		newNZ += nvar - Wcol[Wnz-1];
	}
	Knz = newNZ + Anz + ncon;
	printf("Amount of space required %ld\n", Knz);
	return Knz;
}