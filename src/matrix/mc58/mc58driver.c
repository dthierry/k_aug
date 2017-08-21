
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern void mc58id_(double *CNTL, int *ICNTL);
extern void mc58ad_(int *M, int *N, int* NZ, int *LA, double *A, int *LIRN,
	int *IRN, int *LJCN, int *JCN, double *CNTL, int *ICNTL, int *LIW,
	int *IW, int *INFO, double *RINFO, int *RANK, int *ROWS, int *COLS);

int main(void){
	int m, n, ne, la, lirn, *irn, ljcn, *jcn, liw, *iw, rank, *rows, *cols;
	double *a;
	int info[20], icntl[20];
	double rinfo[10], cntl[10];
	int i, mmn;

	FILE *f_in, *f_out;

	f_in = fopen("jacobi_debug.in", "r");

	if(!f_in){
		printf("I[[MC58DRI]]\tFile not found.\n");
		return 1;
	}
	else{printf("I[[MC58DRI]]\tFile successfully read.\n");}

	fscanf(f_in, "%d\t%d\t%d", &m, &n, &ne);

	if(m < 1 || n < 1 || ne < 1){
		printf("E[[MC58DRI]]\tm n or ne is not admissible.\n");
		return 0;
	}


	printf("Number of rows %d\tNumber of columns %d\t Non-zeroes %d\n", m, n, ne);
	la = 100 * ne;
	lirn = 100 * ne;
	ljcn = 100 * ne;


	printf("I[[MC58DRI]]\tAllocating memory.\n");
	a = (double *)calloc(sizeof(double), la);
	irn = (int *)calloc(sizeof(int), lirn);
	jcn = (int *)calloc(sizeof(int), ljcn);
	
	if(!a || !irn || !jcn){
		printf("E[[MC58DRI]]\tError allocating memory.\n");
	}

	printf("E[[MC58DRI]]\tReading Jacobian matrix..\n");

	for(i=0; i<ne; i++){fscanf(f_in, "%d\t%d\t%lf", irn + i, jcn + i, a + i);}
	fclose(f_in);
	/*for(i=0; i<ne; i++){
		fprintf(stderr, "%d\t%d\t%f\n", *(irn + i), *(jcn + i), *(a + i));
	}*/

	/* LIW ≥ 6*(M+N)+MAX(M,N) */
	mmn = (m > n) ? m:n;
	rows = (int *)calloc(sizeof(int), mmn);
	cols = (int *)calloc(sizeof(int), mmn);
	printf("I[[MC58DRI]]\tMax m n %d\n", mmn);
	
	liw = 6 * (m + n) + mmn + m; /* For good measure */ 
	printf("I[[MC58DRI]]\tCalculated liw %d\n", liw);
	iw = (int *)calloc(sizeof(int), liw);
	icntl[3] = 4;

	mc58id_(cntl, icntl);
	mc58ad_(&m, &n, &ne, &la, a, &lirn,	irn, &ljcn, jcn, cntl, icntl, &liw,
		iw, info, rinfo, &rank, rows, cols);
	for(i = 0; i< 10; i++){
		printf("I[[MC58DRI]]\tINFO value %d i %d tries\n", info[0], i);
		if (info[0] == 0){
			printf("I[[MC58DRI]]\tOperation successful %d\n", info[0]);
			break;
		}

		if(info[0] == -4){
			la = info[4-1];
			lirn = la;
			ljcn = la;

		}
		else if(info[0] == -5){
			lirn = info[5-1];
			la = lirn;
			ljcn = lirn;
		}
		else if(info[0] == -6){
			ljcn = info[6-1];
			la = ljcn;
			lirn = ljcn;
		}

		free(a);
		free(irn);
		free(jcn);
		a = (double *)calloc(sizeof(double), la);
		irn = (int *)calloc(sizeof(int), lirn);
		jcn = (int *)calloc(sizeof(int), ljcn);
		f_in = fopen("jacobi_debug.in", "r");
		for(i=0; i<ne; i++){fscanf(f_in, "%d\t%d\t%lf", irn + i, jcn + i, a + i);}
		fclose(f_in);
		mc58ad_(&m, &n, &ne, &la, a, &lirn,	irn, &ljcn, jcn, cntl, icntl, &liw,
			iw, info, rinfo, &rank, rows, cols);
		}


	f_out = fopen("mc58_out.txt", "w");
	for(i=0; i<rank; i++){
		fprintf(f_out, "%d\t%d\t%d\n", *(rows + i), *(cols + i), rank);
	}
	fclose(f_out);

	free(a);
	free(irn);
	free(jcn);
	free(rows);
	free(cols);
	free(iw);

	return 0;
}


