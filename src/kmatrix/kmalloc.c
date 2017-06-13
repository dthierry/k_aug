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
#include <sys/types.h>
/* Calculates the amount of space for the whole KKT matrix */
fint k_malloc(fint *Wrow, fint *Wcol, real *Wij, fint Wnz, fint nvar, fint ncon,
 fint Anz){

	fint newNZ;
	fint i;
	int j, k;

	
	fint ptr0;
	int miss_nz, miss_row;
	fint Knz;
	ptr0 = 0;
	newNZ = 0;
	miss_row = 0;
	miss_nz = 0;
	
	/* should skip this if Wnz == 0 */
	/* count the number of missing nz */
	/* first loop on vector elements (starts 1) */
	for(i=1; i <= nvar; i++){
		for(j=ptr0; j<Wnz ;j++){
			/* found col */
			if(Wcol[j] == i){
				k = j;
				while(k < Wnz-1 && Wcol[k] == Wcol[j]){
					/* found nz in main diag */
					if(Wrow[k] == i){
						break;
					}
					k++;
				}
				ptr0 = k;
				if(Wrow[k] == i){
					break;
				}
				/* not found nz in main diag */
				else{
					miss_nz++;
					break;
				}
			}
			/* not found col */
			else if (Wcol[j] > i){
				miss_row++;
				break;
			}
			/* else{ 
				printf("Something else happened\n");
			} */
		}
	}

	printf("I[KMATRIX]...\t[KMALLOC]"
			"Number of elements missing in the main diag W%d\n", miss_nz);
	/* printf("Number of elements missing in the main diag W%d\n", miss_nz); */
	
	printf("I[KMATRIX]...\t[KMALLOC]"
			"Number of rows missing W %d\n", miss_row);

	newNZ = miss_row + miss_nz + Wnz;
	
	/* assert(Wcol == NULL); */
	if(Wnz != 0){
		if(Wcol[Wnz-1] < nvar){
			printf("I[KMATRIX]...\t[KMALLOC]"
				"Last element provided is less the square matrix[n_var] for W\n");
			printf("I[KMATRIX]...\t[KMALLOC]"
				"Offset of %ld\n", nvar - Wcol[Wnz-1]);
			newNZ += nvar - Wcol[Wnz-1];
		}
	}
	else{
		printf("I[KMATRIX]...\t[KMALLOC]"
				"No NZ in the hessian detected, set up (n_var) %ld\n", nvar);
		newNZ += nvar;
	}



	Knz = newNZ + Anz + ncon;
	printf("I[KMATRIX]...\t[KMALLOC]Amount of space required %ld\n", Knz);
	return Knz;
}