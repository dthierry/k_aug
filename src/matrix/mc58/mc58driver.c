
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern void mc58id_(double *CNTL, int *ICNTL);
extern void mc58ad_(int *M, int *N, int* NZ, int *LA, double *A, int *LIRN,
	int *IRN, int *LJCN, int *JCN, double *CNTL, int *ICNTL, int *LIW,
	int *IW, int *INFO, double *RINFO, int *RANK, int *ROWS, int *COLS);

int main(void){
	int m, n, ne, la, lirn, *irn, ljcn, *jcn, *iw, rank, *rows, *cols;
	double *a;
      int liw;
	int info[20], icntl[20];
	double rinfo[10], cntl[10];
	int i, j, mmn;

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
	la =   3 * ne;
	lirn = 3 * ne;
	ljcn = 3 * ne;


	printf("I[[MC58DRI]]\tAllocating memory for A matrix with %d.\n", la);
	a = (double *)malloc(sizeof(double)* la);
	irn = (int *)malloc(sizeof(int)* lirn);
	jcn = (int *)malloc(sizeof(int)* ljcn);
      
      memset(a, 0, sizeof(double)*la);
      memset(irn, 0, sizeof(int)*lirn);
      memset(jcn, 0, sizeof(int)*ljcn); 

	
	if(!a || !irn || !jcn){
		printf("E[[MC58DRI]]\tError allocating memory.\n");
	}

	printf("I[[MC58DRI]]\tReading Jacobian matrix..\n");

	for(i=0; i<ne; i++){fscanf(f_in, "%d\t%d\t%lf", irn + i, jcn + i, a + i);}
	fclose(f_in);
      f_in = fopen("debug_.txt", "w");

	for(i=0; i<ne; i++){
	/*	fprintf(stderr, "%d\t%d\t%f\n", *(irn + i), *(jcn + i), *(a + i));*/
		fprintf(f_in  , "%d\t%d\t%f\n", *(irn + i), *(jcn + i), *(a + i));
	}
      fclose(f_in);

      
	/* LIW ≥ 6*(M+N)+MAX(M,N) */
	mmn = (m > n) ? m:n;
	rows = (int *)calloc(sizeof(int), mmn);
	cols = (int *)calloc(sizeof(int), mmn);
	printf("I[[MC58DRI]]\tMax m n %d\n", mmn);
      if (!rows || !cols){
            printf("The memory for rows and cols was not properlly allocated");
            exit(-1);
      }    	
	liw = 6 * (m + n) + mmn; /* For good measure */ 
	printf("I[[MC58DRI]]\tCalculated liw %d\n", liw);
	iw = (int *)calloc(sizeof(int), liw);
	icntl[4-1] = 4;
      cntl[2-1] = 0.5;
      if (!iw){
            printf("The memory was not properly allocated \n");
            exit(-1);
      }

	mc58id_(cntl, icntl);
      cntl[2-1] = 0.8;
	mc58ad_(&m, &n, &ne, &la, a, &lirn,	irn, &ljcn, jcn, cntl, icntl, &liw,
		iw, info, rinfo, &rank, rows, cols);
	for(j = 0; j< 20; j++){
		printf("I[[MC58DRI]]\tINFO value %d i %d tries\n", info[0], j);
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
            
            la = la * 2;
            lirn = lirn * 2;
            ljcn = ljcn * 2;
            
		free(a);
		free(irn);
		free(jcn);
		a = (double *)calloc(sizeof(double), la);
		irn = (int *)calloc(sizeof(int), lirn);
		jcn = (int *)calloc(sizeof(int), ljcn);
		f_in = fopen("jacobi_debug.in", "r");

      	fscanf(f_in, "%d\t%d\t%d", &m, &n, &ne);
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


