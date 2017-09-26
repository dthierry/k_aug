
#include "ma57_driver.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

/*
extern void ma57id_(double *CNTL, int *ICNTL);
extern void ma57ad_(int *N, int *NE, int *IRN, int *JCN, int *LKEEP, int *KEEP,
 int *IWORK, int *ICNTL, int *INFO, double *RINFO);
extern void ma57bd_(int *N, int *NE, double *A, double *FACT, int *LFACT,
 int *IFACT, int *LIFACT, int *LKEEP, int *KEEP, int *IWORK, int *ICNTL,
 double *CNTL, int *INFO, double *RINFO);
extern void ma57cd_(int *JOB, int *N, double *FACT, int *LFACT, int *IFACT, 
	int *LIFACT,	int *NRHS, double *RHS, int *LRHS, double *WORK, int *LWORK,
	int *IWORK, int *ICNTL, int *INFO);
extern void ma57dd_(int *JOB, int *N, int *NE, double *A, int *IRN, int *JCN,
 double *FACT, int *LFACT, int *IFACT, int *LIFACT, double *RHS, double *X,
 double *RESID, double *WORK, int *IWORK, int *ICNTL, double *CNTL,
 int *INFO, double *RINFO);
extern void ma57ed_(int *N, int *IC, int *KEEP, double *FACT, int *LFACT,
 double *NEWFAC, int *LNEW, int	*IFACT, int *LIFACT, int *NEWIFC, int *LINEW,
  int *INFO);*/


/*
int main(void){
	int irn[] = {1, 2};
	int jrn[] = {1, 2};
	double a[] = {1.0, 1.0};
	double b[] = {1.999999, 2.9e-1};
	double x[] = {0, 0};
	ma57_driver(2, 2, irn, jrn, a, 1, b, x);
	return 0;
}

*/

void ma57_driver(int n, int nz, int *i_rn, int *j_rn, double *a, int n_rhs,
 double *rhs, double *x){
	
	double cntl[5];
	int icntl[20];
	int info[40];
	double rinfo[20];

	int i;
	int space_lk; /* for lkeep */
	int *keep;
	int *iwork;
	int lfact, lifact;
	double *fact;
	int *ifact;

	int lfact_new, lifact_new;
	double *fact_new;
 	int *ifact_new;
	int ic;
	int job;

	double *work;
	int lwork;

	int lrhs;

	double rdummy[] = {0.0};
	int idummy[] = {0};

	int _SUC_FACT = 0;

	double const residual_ratio_max = 1e-10;

	int nev;
	double *resid;

	int incx=1; /* for the norm calculation */
	double nrm_x=0, nrm_r=0;

	space_lk = 5 * n + nz + (n > nz ? n:nz) + 42;
	keep = (int *)malloc(sizeof(int) * space_lk);
	memset(keep, 0, sizeof(int)* space_lk);
	iwork = (int *)malloc(sizeof(int) * 5 * n);

	memset(cntl, 0, sizeof(double)*5);
	memset(icntl, 0, sizeof(int)*20);
	
	ma57id_(cntl, icntl);

	cntl[1-1] = 1e-08;


	icntl[5-1] = 2;
	icntl[11-1] = 16;
	icntl[12-1] = 16;

	printf("I[MA57]...\t[Pivot-Selection]"
		"***\n");

	ma57ad_(&n, &nz, i_rn, j_rn, &space_lk, keep, iwork, icntl, info, rinfo);

	if(info[0] != 0){
		printf("I[MA57]...\t[ma57bd_]"
			"ma57ad_: info[0] is not zero; %d \n", info[1-1]);
		exit(-1);
	}
	else{
		printf("I[MA57]...\t[ma57bd_]"
			"ma57ad_: Status Ok.\n");
	}

	printf("I[MA57]...\t[ma57bd_]"
		"ma57ad_: reported value of INFO(9): %d\n",  info[9-1]);
	printf("I[MA57]...\t[ma57bd_]"
		"ma57ad_: reported value of INFO(10): %d\n", info[10-1]);
		
	printf("I[MA57]...\t[Factorize]"
		"***\n");

	lfact = info[9-1] * 2;
	lifact = info[10-1] * 2;
	
	fact = (double *)malloc(sizeof(double) * lfact);
	ifact = (int *)malloc(sizeof(int) * lifact);
	/*emset(fact, 0, sizeof(double) * lfact);*/

	for(i = 0; i < 20; i++){
		ma57bd_(&n, &nz, a, fact, &lfact,	ifact, &lifact, &space_lk, 
			keep, iwork, icntl,	cntl, info, rinfo);
		printf("\n\n");
		if(info[0] != 0){
			printf("W[MA57]...\t[ma57bd_]"
				"ma57bd_: info[0] is not zero; %d \n\n", info[1-1]);
		}
		else{
			printf("I[MA57]...\t[ma57bd_]"
				"ma57bd_: Status Ok.\n\n");
			_SUC_FACT = 1;
		}

		if(_SUC_FACT == 1){break;}
		
		if(info[0] == -3){
			printf("I[MA57]...\t[ma57bd_]"
				"ma57ad_: bad value of INFO(9). Set up new value: %d\n",  info[17-1]);
			ic = 0;
			lfact_new = info[17-1];
			fact_new = (double *)malloc(sizeof(double) * lfact_new);
			memset(fact_new, 0, sizeof(double) * lfact_new);
			assert(fact_new);
			assert(keep);
			/*MA57ED(N,IC,KEEP,FACT,LFACT,NEWFAC,LNEW,IFACT,LIFACT,NEWIFC,LINEW,INFO)*/
			ma57ed_(&n, &ic, keep, fact, &lfact, fact_new, &lfact_new,
				ifact, &lifact, idummy, idummy,	info);
			if(info[0] != 0){exit(-1);}
			free(fact);
			fact = fact_new;
			lfact = lfact_new;
			fact_new = NULL;
		}
		else if(info[0] == -4){
			printf("I[MA57]...\t[ma57bd_]"
				"ma57ad_: bad value of INFO(10). Set up new value: %d\n",  info[18-1]);
			ic = 1;
			lifact_new = info[18-1];
			ifact_new = (int *)malloc(sizeof(int) * lifact_new);
			assert(ifact_new);
			assert(keep);
			ma57ed_(&n, &ic, keep, fact, &lfact, rdummy, idummy,
				ifact, &lifact, ifact_new, &lifact_new,	info);
			if(info[0] != 0){exit(-1);}
			free(ifact);
			ifact = ifact_new;
			lifact = lifact_new;
			ifact_new = NULL;
		}
		else if(info[0] < 0){
			printf("E[MA57]...\t[ma57bd_]"
				"ma57ad_: info(0): %d\n",  info[0]);
			exit(-1);
		}
	}

	nev = info[24-1];
	printf("I[MA57]...\t[Neg_eigenvalues]"
		"%d\n", nev);
	printf("I[MA57]...\t[Solve]"
		"***\n");

	job = 1;

	lwork = n*n_rhs;
	work = (double *)malloc(sizeof(double)*lwork);
	lrhs = n;

	for(i=0; i<n*n_rhs; i++){
		/*Momentarily swap vectors*/
		x[i] = rhs[i];
	}

	ma57cd_(&job, &n, fact, &lfact, ifact, &lifact,	&n_rhs, rhs, &lrhs,
	 work, &lwork,	iwork, icntl, info);
	
	nrm_x = dnrm2_(&n, rhs, &incx);
	printf("I[MA57]...\t[ma57dd_]"
			"ma57bd_:Eucl Norm of solution; %e \n", nrm_x);

	/*printf("I[MA57]...\t[ma57dd_]"
			"ma57bd_: rinfo[9] inf Norm of solution; %e \n", rinfo[9-1]);*/

	if(info[0] != 0){
			printf("W[MA57]...\t[ma57bd_]"
				"ma57bd_: info[0] is not zero; %d \n\n", info[1-1]);
	}

	/*MA57DD(JOB,N,NE,A,IRN,JCN,FACT,LFACT,IFACT,LIFACT,RHS,X,RESID,WORK,IWORK,
	ICNTL,CNTL,INFO,RINFO)*/

	/* he maximum permitted number of steps of iterative refinement */
	icntl[9-1] = 1;
	job = 1;

	resid = (double *)malloc(sizeof(double) * n);

	/*MA57DD(JOB,N,NE,A,IRN,JCN,FACT,LFACT,IFACT,LIFACT,RHS,X,RESID,WORK,IWORK,
	ICNTL,CNTL,INFO,RINFO)*/
	/* It can be done one vector at a time */
	/* swap x and rhs */
	/* try the first vector for now*/

	ma57dd_(&job, &n, &nz, a, i_rn, j_rn, fact, &lfact, ifact, &lifact, x, rhs, resid, work, iwork, 
		icntl, cntl, info, rinfo);

	nrm_r = dnrm2_(&n, resid, &incx);
	printf("I[MA57]...\t[ma57dd_]"
			"ma57bd_:Eucl Norm of the residuals; %e \n", nrm_r);

	/*printf("I[MA57]...\t[ma57dd_]"
				"ma57bd_: info[10] Norm of scaled residuals; %g \n", rinfo[10-1]);*/
	printf("I[MA57]...\t[ma57dd_]"
				"ma57bd_: Ratio of norm of scaled residuals; %g \n", (nrm_r/(nrm_x+nrm_r)));
	if((nrm_r/(nrm_x+nrm_r)) < residual_ratio_max){
			printf("I[MA57]...\t[ma57bd_]"
				"ma57dd_: The norm of residuals is less than max ratio\n");
	}

	for(i=0; i<n*n_rhs; i++){
		x[i] = rhs[i];
	}

	
	
	free(resid);
	free(keep);
	free(iwork);
	free(fact);
	free(ifact);
	free(work);
}
