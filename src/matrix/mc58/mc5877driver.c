
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern void mc58id_(double *CNTL, int *ICNTL);
extern void mc58ad_(int *M, int *N, int* NZ, int *LA, double *A, int *LIRN,
	int *IRN, int *LJCN, int *JCN, double *CNTL, int *ICNTL, int *LIW,
	int *IW, int *INFO, double *RINFO, int *RANK, int *ROWS, int *COLS);

extern void mc77id_(int *icntl, double *cntl);
extern void mc77bd_(int *JOB, int *M, int *N, int *NNZ, int *IRN, int *JCN, double *A, int *IW, int *LIW, double *DW, int *LDW, int *ICNTL, double *CNTL, int *INFO, double *RINFO);


int main(int argc, char *argv[]){
	int m, n, ne, la, lirn, *irn, ljcn, *jcn, *iw, rank, *rows, *cols;
	double *a, *ascaled;
      int *irn0, *jcn0;
      int liw;
	int info[20], icntl[20];
	double rinfo[10], cntl[10];
	int i, icol, irow, j, mmn;

      int liwmc77, ldwmc77, jobmc77;
      int icntlmc77[10], infomc77[10], *iwmc77;
      double cntlmc77[10], *dwmc77, rinfomc77[10];
      

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
	la =   10 * ne;
	lirn = 10 * ne;
	ljcn = 10 * ne;


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

      liwmc77 = mmn * 2;
      ldwmc77 = mmn*n + 2*(m + n);
      
      /* MC77 PART */
      /*Allocate dw, liw*/
      dwmc77  = (double *) malloc(sizeof(double)*ldwmc77);    
      iwmc77 =  (int *)malloc(sizeof(int)*liwmc77);  

      if(!dwmc77 || !iwmc77){
            printf("E[[MC5877DRI]]\tError during allocation of dw and iw (MC77).\n");
            exit(-1);
      }
      memset(dwmc77, 0, sizeof(double)*ldwmc77);
      memset(iwmc77, 0, sizeof(int)   *liwmc77);

      mc77id_(icntlmc77, cntlmc77);

      jobmc77 = 0;
      mc77bd_(&jobmc77, &m, &n, &ne, irn, jcn, a, iwmc77, &liwmc77, dwmc77, &ldwmc77, icntlmc77, cntlmc77, infomc77, rinfomc77);  

      printf("I[[MC5877DRI]]\tICNTL_MC77(7) = %d\n", icntlmc77[7-1]);
      printf("I[[MC5877DRI]]\tICNTL_MC77(6) = %d\n", icntlmc77[6-1]);
      printf("I[[MC5877DRI]]\tRINFO_MC77(1) = %f\n", rinfomc77[1-1]);
      printf("I[[MC5877DRI]]\tRINFO_MC77(2) = %f\n", rinfomc77[2-1]);

      f_out = fopen("mc77out.txt", "w");

      for(i=0; i < m; i++) {fprintf(f_out, "%d\t%f\n", i, dwmc77[i]);}
      fprintf(f_out,"\n\n");
      if(icntlmc77[6-1] == 0){
            for(i=0; i < n; i++) {fprintf(f_out, "%d\t%f\n", i, dwmc77[m + i]);}
      }
      fclose(f_out); 
      
      if(argc > 1){
            printf("I[[MC5877]]\tScaling deactivated\n");
      }
      else{
            printf("I[[MC5877]]\tScaling activated\n");

      /* MATRIX SCALING */
      for(i=0; i< ne; i++){
            irow = irn[i];
            icol = jcn[i];
            a[i] = a[i]*(1/dwmc77[icol-1]);
      }
      if(icntlmc77[6-1] == 0){
            for(i=0; i< ne; i++){
                  irow = irn[i];
                  a[i] = a[i]*(1/dwmc77[m - 1 + irow]);
            }
      }
      else{ 
            for(i=0; i< ne; i++){
                  irow = irn[i];
                  icol = jcn[i];
                  a[i] = a[i]*(1/dwmc77[irow-1]);
            }
      }
      }
      free(dwmc77);
      free(iwmc77);

      ascaled = (double *)calloc(sizeof(double), ne);
      irn0 = (int *)calloc(sizeof(int), ne);
      jcn0 = (int *)calloc(sizeof(int), ne);



      for(i=0; i<ne; i++){
            ascaled[i] = a[i];
            irn0[i] = irn[i];
            jcn0[i] = jcn[i];
      }

      f_in = fopen("debug_scaled0.txt", "w");
	for(i=0; i<ne; i++){
	/*	fprintf(stderr, "%d\t%d\t%f\n", *(irn + i), *(jcn + i), *(a + i));*/
		fprintf(f_in  , "%d\t%d\t%f\n", *(irn + i), *(jcn + i), *(a + i));
	}
      fclose(f_in);

	rows = (int *)calloc(sizeof(int), mmn);
	cols = (int *)calloc(sizeof(int), mmn);
	printf("I[[MC5877DRI]]\tMax m n %d\n", mmn);
      if (!rows || !cols){
            printf("The memory for rows and cols was not properlly allocated");
            exit(-1);
      }    	
	liw = 6 * (m + n) + mmn; /* For good measure */ 
	printf("I[[MC5877DRI]]\tCalculated liw %d\n", liw);
	iw = (int *)calloc(sizeof(int), liw);
      if (!iw){
            printf("The memory was not properly allocated \n");
            exit(-1);
      }

	mc58id_(cntl, icntl);
      cntl[2-1] = 0.9;
      icntl[4-1] = 3;
      icntl[8-1] = 32;
	mc58ad_(&m, &n, &ne, &la, a, &lirn,	irn, &ljcn, jcn, cntl, icntl, &liw,
		iw, info, rinfo, &rank, rows, cols);
	for(j = 0; j< 1000; j++){
		printf("I[[MC5877DRI]]\tINFO value %d i %d tries\n", info[0], j);
		if (info[0] == 0){
			printf("I[[MC5877DRI]]\tOperation successful %d\n", info[0]);
			break;
		}

		if(info[0] == -4){
			la = info[4-1];
/*			lirn = la;
			ljcn = la;*/

		}
		else if(info[0] == -5){
			lirn = info[5-1];
/*			la = lirn;
			ljcn = lirn;*/
		}
		else if(info[0] == -6){
			ljcn = info[6-1];
/*			la = ljcn;
			lirn = ljcn;*/
		}
            else{exit(-1);}
            
            
            la = la*2 ;
            lirn = lirn*2 ;
            ljcn = ljcn*2 ;
          	mc58id_(cntl, icntl);
            cntl[2-1] = 0.9;
            icntl[4-1] = 3;

		free(a);
		free(irn);
		free(jcn);
		a = (double *)malloc(sizeof(double)*la);
		irn = (int *)malloc(sizeof(int)*lirn);
		jcn = (int *)malloc(sizeof(int)*ljcn);

 /*         memset(a, 0, sizeof(double)*la);
            memset(irn, 0, sizeof(int)*lirn);
            memset(jcn, 0, sizeof(int)*ljcn);
  **/

		/*for(i=0; i<ne; i++){
                  a[i] = ascaled[i];
                  irn[i] = irn0[i];
                  jcn[i] = jcn0[i];
            }*/

      f_in = fopen("debug_scaled0.txt", "r");
      f_out = fopen("debug_scaled.txt", "w");
	for(i=0; i<ne; i++){
		fscanf(f_in  , "%d\t%d\t%lf\n", (irn + i), (jcn + i), (a + i));
		fprintf(f_out, "%d\t%d\t%f\n", *(irn + i), *(jcn + i), *(a + i));

	}
      fclose(f_in);
      fclose(f_out);


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
      free(ascaled);
      free(irn0);
      free(jcn0);

	return 0;
}


