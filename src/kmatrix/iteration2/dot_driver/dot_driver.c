/* @source untitled.c
** beta 01
** Month dayth, 20yy
** @author: David Thierry (dmolinat@andrew.cmu) dav0@lb2016-1
	READ THE DATA FOR S HAT
	NEED SUFFIX FOR DELTA P

********************************************************************************

@fun_name ********************************************
**
** Description
** Description
**
** @param [r] 
** @param [r] 
** @param [r] 
** @param [r] 
** @param [r] 
** @param [r] 
** @param [r] 
** @return something
*******************************************************************************/

#include "asl.h"
#include "getstub.h"
#include <stdio.h>
#include <stdlib.h>
#include "sens_update_driver_dot.h"

static Option_Info Oinfo;
static keyword *keywds;

int main(int argc, char **argv){
	ASL *asl;
	unsigned i, j, k;
	unsigned nvar, ncon;
	unsigned _n_dof;
	unsigned n_srow;

	char *s=NULL;
	FILE *f;
	FILE *rh_txt;
	FILE *f_out;

	double *x=NULL, *lambda=NULL;
	double *s_hat_T=NULL;
	double *s_hat_=NULL;

	double *npdp=NULL;
	double *u_star=NULL;
	double *u_star0=NULL;

	int *u_arr=NULL;
	
	SufDecl *suf_ptr=NULL;

	char _suf1[] = {"dof_v"};
	char _suf2[] = {"npdp"};
	SufDesc *suf1=NULL;
	SufDesc *suf2=NULL;

  char ter_msg[] = {"I[[SENS_DOT]]...\t[MAIN]"
	"All done it seems."};

	char _sname[] = {"[[SENS_DOT]]"};

	Oinfo.sname = _sname;
	Oinfo.bsname = NULL;
	Oinfo.opname = NULL;
	Oinfo.keywds = NULL;
	Oinfo.n_keywds = nkeywds;
	Oinfo.flags = 0;
	Oinfo.version = NULL;
	Oinfo.usage = NULL;
	Oinfo.kwf = NULL;
	Oinfo.feq = NULL;
	Oinfo.n_options = 0;
	Oinfo.driver_date = 0;
	Oinfo.wantsol = 0;
	Oinfo.nS = 0;
	Oinfo.S = NULL;
	Oinfo.uinfo = NULL;
	Oinfo.asl = NULL;
	Oinfo.eqsign = NULL;
	Oinfo.n_badopts = 0;
	Oinfo.option_echo = 0;
	Oinfo.nnl = 0;

	asl = ASL_alloc(ASL_read_fg);
	s = getstops(argv, &Oinfo);

	if(!s){
		printf("I[[SENS_DOT]]...\t[MAIN]"
			"No input.\n");
		ASL_free(&asl);
		return 0;
	}
	else{
		printf("I[[SENS_DOT]]...\t[MAIN]"
			"File read succesful.\n");
	}

	/* Two suffixes */
	suf_ptr = (SufDecl*)M1alloc(sizeof(SufDecl) * 2);

	
	/* consider using other suffix instead */
	suf_ptr->name = _suf1;
  suf_ptr->table = 0;
  suf_ptr->kind = ASL_Sufkind_var;
  suf_ptr->nextra = 0;

  (suf_ptr+1)->name = _suf2;
  (suf_ptr+1)->table = 0;
  (suf_ptr+1)->kind = ASL_Sufkind_con | ASL_Sufkind_real;
  (suf_ptr+1)->nextra = 0;

  suf_declare(suf_ptr, 2);

	f = jac0dim(s, (int)strlen(s));
	nvar = (unsigned)n_var;
	ncon = (unsigned)n_con;
	n_srow = nvar + ncon;
	printf("I[[SENS_DOT]]...\t[MAIN]"
		"Number of variables %d\n", (int)nvar);
	printf("I[[SENS_DOT]]...\t[MAIN]"
		"Number of constraints %d\n", (int)ncon);
	printf("I[[SENS_DOT]]...\t[MAIN]"
		"Number of rows (primal-dual) %d\n", (int)n_srow);

	/* Primal and dual */
	X0 = M1alloc(sizeof(double) * nvar);
	pi0= M1alloc(sizeof(double) * ncon);

	fg_read(f, 0);

	x = X0;
	lambda = pi0;

	suf1 = suf_get(_suf1, ASL_Sufkind_var); /* rh_name */
	suf2 = suf_get(_suf2, ASL_Sufkind_con); /* npdp */

	if(!(suf1->u.i)){
		printf("E[[SENS_DOT]]...\t[MAIN]"
			"No \"%s\" suffix declared. Exiting.\n\n", _suf1);
		ASL_free(&asl);
		exit(-1);
	}
	if(!(suf2->u.r)){
		printf("E[[SENS_DOT]]...\t[MAIN]"
			"No \"%s\" suffix declared. Exiting.\n\n", _suf2);
		ASL_free(&asl);
		exit(-1);
	}
	npdp = (double *)malloc(sizeof(double)*n_srow);
	memset(npdp, 0, sizeof(double)*n_srow);

	for(i=0; i<ncon; i++){
		*(npdp + nvar + i) = suf2->u.r[i];
	}
	
	u_arr = (int *)malloc(sizeof(int)*nvar);

	_n_dof = 0; /* count_dof */
	for(i=0; i<nvar; i++){
		if(*((suf1->u.i) + i) != 0){
			u_arr[_n_dof] = (int)i;
			_n_dof++;
		}
	}

	u_star = (double *)malloc(sizeof(double) * _n_dof);
	u_star0= (double *)malloc(sizeof(double) * _n_dof);
	for(i=0; i<_n_dof; i++){
		u_star[i] = X0[u_arr[i]];
	}

	printf("I[[SENS_DOT]]...\t[MAIN]"
		"Number of dof detected  %u\n",_n_dof);	

	rh_txt = fopen("../dot_in.in", "r");

	

	if(!rh_txt){
		printf("I[[SENS_DOT]]...\t[MAIN]File not found.\n");
		ASL_free(&asl);
		free(u_arr);
		free(u_star);
		exit(-1);
	}

	s_hat_T = (double *)malloc(sizeof(double) * n_srow * _n_dof);
	s_hat_  = (double *)malloc(sizeof(double) * n_srow * _n_dof);
	

	for(i=0; i<_n_dof; i++){
		for(j=0; j<n_srow; j++){
			fscanf(rh_txt, "%lf", (s_hat_T + n_srow * i + j));
		}
	}


	for(i=0; i<n_srow; i++){
		for(j=0; j<_n_dof; j++){
			s_hat_[i*_n_dof + j] = s_hat_T[j*n_srow + i];
		}
	}

	
	fclose(rh_txt);
	
	memcpy(u_star0, u_star, sizeof(double)*_n_dof); /* create a copy of u for debugging*/

	sens_update_driver_dot((fint)n_srow, (fint)_n_dof, s_hat_T, npdp, u_star);
	
	f_out = fopen("dot_out.out", "w");

	for(i=0; i<_n_dof;i++){
		fprintf(f_out, "%.g\t%.g\t%.g\n", u_star0[i], u_star[i], u_star0[i] - u_star[i]);
	}


	for(i=0; i<_n_dof; i++){
		X0[u_arr[i]] = u_star[i]; /* update primal */
	}

	write_sol(ter_msg, X0, pi0, 0);
	printf("I[[SENS_DOT]]...\t[MAIN]"
		"Done.\n");
	fclose(f_out);

	free(npdp);
	free(s_hat_T);
	free(s_hat_);
	free(u_arr);
	free(u_star);
	ASL_free(&asl);
	return 0;

}