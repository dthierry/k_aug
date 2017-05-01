#include "csort.h"

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
    // start the array at row 1
    j = 1;
    k = 0;
    // create array that contains copy of pointers
    for(i=0; i<nz; i++){
    	ttemp[i].c = col+i;
    	ttemp[i].a = aij+i;
    }
    // assign helper pointer to the array
    t = ttemp;

    // do the actual reorder
    for(i=0; i<nz; i++){
    	if (row[i] > j){
            if(row[i] > j + 1){
                fwarn = 1;
                printf("Warning missing row in A by: %ldrows\n",(row[i]-j));
            }
    		qsort(t, k, sizeof(temp), compf);
            /*
            printf("k value %d\n", k);
            for(int l=0; l<k; l++){
                printf("ordered %ld\t", (t+l)->c[0]);
                printf("ordered %f\n", (t+l)->a[0]);
            }
    		printf("\n");
            */
            t = t + k;
    		k = 0;
    		j++;
        }
    	k++;
    }
    
 	qsort(t, k, sizeof(temp), compf);
    /*
     for(int l=0; l<k; l++){
        printf("ordered %ld\t", (t+l)->c[0]);
        printf("ordered %f\n", (t+l)->a[0]);
      }
    */
    // move helper pointer to last position in the array
    t = t + k-1;
 
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

