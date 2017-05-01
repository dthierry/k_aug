#include <stdio.h>
#include <stdlib.h>

#define fint long
#define real double

typedef struct temp{
	fint *c;
	real *a;
}temp;

int compf(const void *t1, const void *t2);
int ccolum(fint *row, fint *col, real *aij, fint nz, fint nrow);



int main(){
 fint cols[] ={439,33,1355,440,34,1356,441,35,1357};
 real aij[] = {-1,202466.143011526,11105.8179573795,-1,203932.427131064,10832.008939236,-1,204850.757715235,10671.5147167406};
 fint rows[] = {1, 1, 1, 2, 3, 4, 5, 5, 6};
 fint nrow = 3;
 fint nz = 9;
 ccolum(rows, cols, aij, nz, nrow);
}

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
 

    return 0;
}

