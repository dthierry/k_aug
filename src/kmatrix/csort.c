#include "csort.h"

int compf(const void *t1, const void *t2)
{
	temp firstt = *(const temp *)t1;
	temp secondt = *(const temp *)t2;
	int fst = (int) (firstt.c[0]);
	int scd = (int) (secondt.c[0]);
	return (fst - scd);
}



int ccolum(fint *row, fint *col, real *aij, fint nz, fint nrow){
	

    temp ttemp[nz];
    temp *t;
    fint cbak[nz];
    real abak[nz];
    
    int i, j, k;

    j = 1;
    k = 0;
    
    for(i=0; i<nz; i++){
    	ttemp[i].c = col+i;
    	ttemp[i].a = aij+i;
    }
    t = ttemp;
    for(i=0; i<nz; i++){
    	if (row[i] > j){
    		qsort(t, k, sizeof(temp), compf);
        printf("k value %d\n", k);
          for(int l=0; l<k; l++){
            printf("ordered %d\t", (t+l)->c[0]);
            printf("ordered %f\n", (t+l)->a[0]);

          }
    		printf("\n");
        t = t + k;

    		k = 0;
    		j++;
    		
    	}
    	k++;

    }
 		qsort(t, k, sizeof(temp), compf);
      for(int l=0; l<k; l++){
            printf("ordered %d\t", (t+l)->c[0]);
            printf("ordered %f\n", (t+l)->a[0]);
      }
 
    // set pointers to the desired value
     //cbak = col;
     //abak = aij;
     for(i=0; i<nz; i++){
    	cbak[i] = (t+i)->c[0];
    	abak[i] = (t+i)->a[0];
      }
    for(i=0; i<nz; i++){
    col[i] = cbak[i];
    aij[i] = abak[i];
    }
    
       
    return 0;
}

