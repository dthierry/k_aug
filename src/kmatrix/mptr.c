#include <stdio.h>
#include <stdlib.h>

extern void matrixptr_(int *a, int *n);

int main(){
	int *a;
	a = (int *)malloc(sizeof(int)*10);
	for(int i=0; i<10; i++){
		a[i] = 4;
		if(i == 2){
			a[i] = 12;
		}
	}
	int n = 10;
	matrixptr_(a, &n);

	for(int i=0; i<10; i++){
		printf("a[i] = %d\n", a[i]);
	}

	free(a);
	return 0;
}