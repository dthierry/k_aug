#include "csort.h"
#include <stdio.h>
int compf(const void *t1, const void *t2)
{
	temp firstt = *(const temp *)t1;
	temp secondt = *(const temp *)t2;
	int fst = (int) (firstt.c[0]);
	int scd = (int) (secondt.c[0]);
	return (fst - scd);
}



int sortcol(fint *row, fint *col, real *aij, fint nz, fint nrow){
	

    temp ttemp[nz];
    temp *t;
    fint cbak[nz];
    real abak[nz];
    unsigned char fwarn=0;
    int i, j, k;
    // start the array at row[0]
    j = 0;
    k = 1;

    // create array that contains copy of pointers
    for(i=0; i<nz; i++){
    	ttemp[i].c = col+i;
    	ttemp[i].a = aij+i;
    }
    // assign helper pointer to the array
    t = ttemp;

    // do the actual reorder
    for(i=0; i < (nz-1); i++){
    	if (row[i] < row[i+1]){
            // offset greater than 1
            if(row[i] - row[i+1] >  1){
                fwarn = 1;
                printf("Warning missing row in A by: %drows\n",(row[i]-row[i+1]));
            }
            // if k == 1 there is no reordering to be done
            if(k > 1){qsort(t, k, sizeof(temp), compf);}
            t += k;
    		k = 1;
    		j++;
        }
        else{
            k++;
        }
    }
    // sort the last row
 	if(k > 1){qsort(t, k, sizeof(temp), compf);}

    // move helper pointer to last position in the array
    t = t + k-1;
    printf("Rows in the matrix: %d\n", nrow);
    printf("Rows processed: %d\n plus 1", j);

    // set pointers to the desired value   
    // move from last to first
    for(i=0; i<nz; i++){
    	cbak[nz-1 - i] = (t-i)->c[0];
    	abak[nz-1 - i] = (t-i)->a[0];
        //printf("column %ld\n", cbak[nz - i]);
      }
    // update back the values of col and aij
    for(i=0; i<nz; i++){
    col[i] = cbak[i];
    aij[i] = abak[i]; 
    }
    
       
    return 0;
}

