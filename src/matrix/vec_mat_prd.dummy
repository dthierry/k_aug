#include <stdlib.h>
#include <stdio.h>

/*
 dgemv(character TRANS, integer M, integer N, double precision ALPHA, double precision, dimension(lda,*) 	A, integer LDA, double precision, dimension(*) 	X)
 integer INCX, double precision BETA, double precision, dimension(*) Y,integer INCY
*/
void dgemv_(char *TRANS, int *M, int *N, double *ALPHA, double *A, int *LDA, double *X, int *INCX, double *BETA, double *Y, int *INCY);

int main(void){
	double *A = NULL;
	double *x = NULL;
	int i, n;
	FILE *ptr_file;

	char t = 'N';
	int M = 10;
	int N = 2;
	double ALPHA = 1.0;
	int LDA = 10;
	int INCX = 1;
	double BETA = 1.0;
	double *Y = NULL;
	int INCY = 1;

	n = 2;
	x = (double *)calloc(N, sizeof(double));
	Y = (double *)calloc(M, sizeof(double));
	A = (double *)calloc(N * M, sizeof(double));

	for(i=0; i<M; i++){
		Y[i] = 0.0;
	}

	ptr_file = fopen("x.txt", "r");

	for(i=0; i<n; i++){
		fscanf(ptr_file, "%lf", (x + i));
	}

	fclose(ptr_file);

	for(i=0; i<n; i++){
		printf("x %f\n", x[i]);
	}
	
	ptr_file = fopen("a.in", "r");

	for(i=0; i<M; i++){
		/* fscanf(ptr_file, "%lf %lf", (A + i*n), (A + i*n + 1));	 */
		fscanf(ptr_file, "%lf %lf", (A + i), (A + M + 1));
	}

	fclose(ptr_file);

	for(i=0; i<M; i++){
		printf("a %f %f\n", A[i], A[M + 1]);
	}

	dgemv_(&t, &M, &N, &ALPHA, A, &LDA, x, &INCX, &BETA, Y, &INCY);	

	for(i=0; i<M; i++){
		printf("y %f\n", Y[i]);
	}

	free(x);
	free(Y);
	free(A);

	return 0;
}