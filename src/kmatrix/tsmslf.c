#include <stdio.h>


int sumefun(int x);


void main(){
	int i;
	printf("Enter value: ");
	scanf("%d", &i);
	sumefun(i);
}



int sumefun(int x){
	int k[x];
	for(int i=0; i<x; i++){
		k[i] = 0;
		printf("The value of k %d\n", k[i]);
	}
}