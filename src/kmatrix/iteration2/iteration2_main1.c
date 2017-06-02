/* @source kmatrix_b0.c
** beta 0
** April 18th, 2017
** @author: David Thierry (dmolinat@andrew.cmu) dav0@lb2016-1

********************************************************************************

@main ********************************************
**
** Reads nl file, allocates data structures, calls assembling funcs
**
** @param [r] stub
** @param [r] KWs
** @@
*******************************************************************************/

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "getstub.h"
#include "w_append_nz.h"
#include "k_assemble_cc.h"
#include "kmalloc.h"
#include "mc30_driver.h"
#include "pardiso_driver.h"

#include "get_jac_asl_aug.h"
#include "get_hess_asl_aug.h"

#include "find_inequalities.h"
#include "assemble_corrector_rhs.h"
#include "assemble_rhsd_red_hess.h"
/*#include "sens_update_driver.h"*/
#include "suffix_decl_hand.h"
#include "csr_driver.h"

static int dumm = 1;
static I_Known dumm_kw = {2, &dumm};
static int n_rhs = 0;
static int n_dof = 0;
static int l_over = 0;
static I_Known l_over_kw = {1, &l_over};
static char name1[] = {"smth"};
static char _n_dofopt_[] = {"n_dof"};
static char _n_rhsopt_[] = {"n_rhs"};
static char _no_lambdaopt_[] = {"no_lambda"};

/*static char dof_v[] = {"dof_v"};*/

/* keywords */
static keyword keywds[] = {
  KW(name1 , IK_val, &dumm_kw, name1),
  KW(_n_dofopt_ , I_val, &n_dof, _n_dofopt_),
  KW(_n_rhsopt_ , I_val, &n_rhs, _n_rhsopt_),	
  KW(_no_lambdaopt_ , IK_val, &l_over_kw, _no_lambdaopt_)
};
static char banner[] = {"[KMATRIX] written by DT\n\n"};
static char _k_[] = {"K_augmented"};
static char _k_o_[] = {"K_augmented options"};
static Option_Info Oinfo = 
{_k_, banner, _k_o_, keywds, nkeywds};
/*Remember to set reference back*/


int main(int argc, char **argv){
	ASL *asl;
	FILE *f;
	FILE *f_debug;
	/* SufDesc *some_suffix; */
	int i, j;
	int nnzw; /* let's try this */
	real *x=NULL, *lambda=NULL;
	char *s=NULL;
	SufDesc *var_f=NULL;
	SufDesc *con_f=NULL;
	SufDesc **rhs_ptr=NULL;
	SufDecl *suf_ptr=NULL;
	SufDecl *some_suffix=NULL;

	fint *Acol=NULL, *Arow=NULL;
	real *Aij=NULL;
	fint *Wcol=NULL, *Wrow=NULL;
	real *Wij=NULL;
	fint *Kcol=NULL, *Krow=NULL;
	real *Kij=NULL;
	fint *Kr_strt=NULL;
	fint K_nrows;

	real *S_scale=NULL;
	fint nzK;

	fint *Wc_t=NULL, *Wr_t=NULL;
	real *Wi_t=NULL;
	fint wa_tsp;

	real *x_=NULL;
	real *dp_=NULL;
	/*real *gf= NULL; */
	real *s_star = NULL;
	int *c_flag=NULL;
	char **rhs_name=NULL;
	char **reg_suffix_name=NULL;
	int *nz_row_w=NULL;
	int *nz_row_a=NULL;
	int *md_off_w=NULL;
	int miss_nz_w;

	/* the pointer to array that will contain */
	/* the rhs */
	real *rhs_baksolve = NULL;

  FILE *somefile;
  fint nerror;
  int nzW, nzA;
  int *hr_point = NULL;
  int *positions_rh = NULL;

  char ter_msg[] = {"I[KMATRIX]...\t[KMATRIX_ASL]"
	"All done it seems."};

	char _dof_s[] = {"dof_v"};
	char _rh_name_s[] = {"rh_name"};


	/* The memory allocation for asl data structure */
	asl = ASL_alloc(ASL_read_pfgh);

	s = getstops(argv, &Oinfo);

	if (!s) {
		printf("W[KMATRIX]...\t[KMATRIX_ASL]"
			"No input\n");
		return 1;
	}
	else {
		printf("I[KMATRIX]...\t[KMATRIX_ASL]"
			"File read succesfull\n");
				}

	if (n_rhs == 0){
		printf("W[KMATRIX]...\t[KMATRIX_ASL]"
			"No n_rhs declared\n");
	}
	if (n_dof == 0){
		printf("W[KMATRIX]...\t[KMATRIX_ASL]"
			"No n_dof declared\n");
	}
	if (n_rhs == 0 && n_dof == 0){
		printf("W[KMATRIX]...\t[KMATRIX_ASL]"
			"Neither n_rhs nor n_dof specified... exiting\n");
		return 1;
	}

	reg_suffix_name =   (char **)malloc(sizeof(char *) * 2);
	reg_suffix_name[0] = (char *)malloc(sizeof(char) * 16 );
	reg_suffix_name[1] = (char *)malloc(sizeof(char) * 16 );
	reg_suffix_name[0][0] = '\0';
	reg_suffix_name[1][0] = '\0';
	strcat(reg_suffix_name[0], _dof_s);
	strcat(reg_suffix_name[1], _rh_name_s);


	if(n_rhs > 0){
		suf_ptr = (SufDecl *)malloc(sizeof(SufDecl)*(n_rhs+2));
		rhs_name = (char **)malloc(sizeof(char *)*n_rhs);
		for(i=0; i<n_rhs; i++){
			rhs_name[i] = (char *)malloc(sizeof(char) * 32); /* 32 bit long digit;
			 why not?*/
		}
		suffix_decl_hand(n_rhs, suf_ptr, reg_suffix_name, rhs_name);
	}
	else{
		suf_ptr = (SufDecl *)malloc(sizeof(SufDecl)*2);
		suffix_decl_hand(n_rhs, suf_ptr, reg_suffix_name, rhs_name);
	}


	printf("I[KMATRIX]...\t[KMATRIX_ASL]"
	"Number of Right hand sides %d\n", n_rhs);
	
	
	for(i=0; i<n_rhs; i++){
		printf("NAME %s\n", rhs_name[i]);
	}

	/* Declare suffixes */
	if(n_rhs > 0){suf_declare(suf_ptr, (n_rhs + 2));}
	else{suf_declare(suf_ptr, (2));}
	

	/* dhis bit setups ASL components e.g. n_var, n_con, etc. */
	f = jac0dim(s, (fint)strlen(s));


	printf("I[KMATRIX]...\t[KMATRIX_ASL]"
		"Number of Right hand sides %d\n", n_rhs);
	printf("I[KMATRIX]...\t[KMATRIX_ASL]"
		"Number of variables %d\n", n_var);
	printf("I[KMATRIX]...\t[KMATRIX_ASL]"
		"Number of constraints %d\n", n_con);

	if (n_dof > (n_var - n_con)){
		printf("E[KMATRIX]...\t[KMATRIX_ASL]"
			"No valid number of n_dof declared\n");
		return -1;
	}
	else if (n_dof < (n_var - n_con) ){
		printf("I[KMATRIX]...\t[KMATRIX_ASL]"
			"n_dof is less than n_var - n_con \n");
	}
	else if (n_dof == (n_var - n_con) ){
		printf("I[KMATRIX]...\t[KMATRIX_ASL]"
			"n_dof exactly n_var - n_con \n");
	}

	x = X0 = M1alloc(n_var*sizeof(real));

	/* Want initial guess lambda */
	want_xpi0 = 2;
  /* I think this allocates _pi0 if = 2 although its behavior is still unknown*/
	/* need to do part of changing sign for y */

	pfgh_read(f, ASL_findgroups);
	/*printf("%d LBC %f\n",0, Urhsx[2*0]); */
	/*find_bounds(n_con, LUv); */
	
	/* change the name of the multipliers to lambda */
	lambda = pi0;

	if(lambda==NULL && l_over == 0){
		printf("E[KMATRIX]...\t[KMATRIX_ASL]"
	"Constraint Multipliers not declared(suffix dual), abort\n");
		for(i=0; i<n_rhs; i++){
			free(rhs_name[i]);
		}
		free(rhs_name);
		ASL_free(&asl);
		free(suf_ptr);
		exit(-1);
	}


	/*con_f = suf_get("con_flag", ASL_Sufkind_con);*/
	var_f = suf_get(reg_suffix_name[0], ASL_Sufkind_var);

	if(var_f->u.r == NULL && var_f->u.i == NULL){
    fprintf(stderr, "E[KMATRIX]...\t[KMATRIX_ASL]"
    	"u.i empty, no var_flag declared.\n");
    exit(-1);
	}
	/*
	if(con_f->u.i == NULL){
		fprintf(stderr, "E[KMATRIX]...\t[KMATRIX_ASL]"
			"u.i empty, no con_flag declared.\n");
		return -1; 
	}
	*/

	rhs_ptr = (SufDesc **)malloc(sizeof(SufDesc *) * n_rhs);

	for(i=0; i < n_rhs; i++){
   *(rhs_ptr + i)= suf_get(rhs_name[i], ASL_Sufkind_con|ASL_Sufkind_real);
  	if((*(rhs_ptr + i))->u.r == NULL){
		  fprintf(stderr, "E[KMATRIX]...\t[KMATRIX_ASL]"
		  	"u.r empty, no rhs values declared for rhs_%d.\n", i);
		  exit(-1);
  	}
	}
	
	c_flag = (int *)malloc(sizeof(int) * n_con); /* Flags for ineq or equalities*/
	
	/*constraints flags */
	find_ineq_con(n_con, LUrhs, c_flag);

	/* Row and colum for the triplet format A matrix */
	/* size of the number of nonzeroes in the constraint jacobian */
	Acol = (fint *)malloc(sizeof(fint)*nzc);
	Arow = (fint *)malloc(sizeof(fint)*nzc);
	Aij  = (real *)malloc(sizeof(real)*nzc);

	/* TO DO:
	CHECH WHICH VARIABLES HAVE BOUNDS
	LOOK FOR BOUND MULTIPLIERS
	CHECH WHICH VARIABLES HAVE MULTIPLIERS
	ASSEMBLE SIGMA
	ADD SIGMA TO HESSIAN

	assemble csr or coordinate
	 */ 
	nerror = 0;
	printf("I[KMATRIX]...\t[KMATRIX_ASL]"
		"Nonzeroes in the sparse Jacobian %d\n", nzc);

	get_jac_asl_aug (asl, x, Acol, Arow, Aij, n_var, n_con, nzc, &nerror, &nz_row_a);
	get_hess_asl_aug(asl, x, &Wcol, &Wrow, &Wij, n_var, n_con, n_obj,
	 &nnzw, lambda, &nerror, &nz_row_w, &md_off_w, &miss_nz_w);
	assert(nerror == 0);
	nzA = nzc; 
	nzW = nnzw + miss_nz_w; /* Recomputed number of nz in the Hessian-of-theLagrange */


	csr_driver((int)n_var, (int)n_con, nzW, nzA, nz_row_w, nz_row_a,
		(int*)Wrow, (int*)Wcol, Wij, (int*)Arow, (int*)Acol, Aij, 
		&Krow, &Kcol, &Kij, &Kr_strt);

	K_nrows = n_var + n_con; /* Number of rows of the KKT matrix (no ineq) */
	nzK = nzA + nzW + n_con; /* NZ in the kkt matrix (for pardiso, yes)*/
	assert(Krow != NULL);

	S_scale = (real *)calloc(sizeof(real), K_nrows);
	mc30driver(K_nrows, nzK, Kij, Krow, Kcol, S_scale);

	/* */
  /* by [rhs][n] */
  rhs_baksolve = (real *)calloc(K_nrows * (n_dof), sizeof(real));
  x_           = (real *)calloc(K_nrows * (n_dof), sizeof(real));
  hr_point		 = (int *)malloc(n_dof * sizeof(int));
  positions_rh = (int *)malloc(n_var * sizeof(int));

  s_star			 = (real *)calloc(K_nrows, sizeof(real));
   /* array that contains the primal dual update */

  /* Primal-dual vector */
  for(i=0; i<K_nrows; i++){
  	if(i<n_var){
  		s_star[i] = x[i];
  	}
  	else{
  		s_star[i] = lambda[i - n_var];
  	}
  }
  
  somefile = fopen("primal_dual.txt", "w");
  
  for(i=0; i<K_nrows; i++){
    fprintf(somefile, "\t%f\n", s_star[i]);
  }
  fclose(somefile);


	/* */
	/*assemble_rhsds(n_rhs, K_nrows, rhs_baksolve, dp_, n_var, n_con, rhs_ptr); */
  assemble_rhs_red_hess(rhs_baksolve, n_var, n_con, n_dof, var_f, hr_point);
  
	
  /* scale matrix */
  for(i=0; i< nzK; i++){
  	Kij[i] = Kij[i] * exp(S_scale[Kcol[i]-1]) * exp(S_scale[Krow[i]-1]);
  }
  
  
  /* v3 */
  for(i=0; i< n_dof; i++){
  	for(j=0; j < K_nrows; j++){
  		*(rhs_baksolve + i*K_nrows + j) = 
  		*(rhs_baksolve + i*K_nrows + j) * exp(S_scale[j]);
  	}
  }

	somefile = fopen("rhs_sens_scaled", "w");


 	/* v3 */
 	for(j=0; j < K_nrows; j++){
 		for(i=0; i < n_dof; i++){
 			fprintf(somefile, "\t%f\t", *(rhs_baksolve + i*K_nrows + j) );
 		}
 		fprintf(somefile, "\n");
 	}
 	fclose(somefile);

 
  /* factorize the matrix */
	pardiso_driver(Kr_strt, Kcol, Kij, K_nrows, nzK, n_dof, rhs_baksolve, x_);
      
  printf("I[KMATRIX]...\t[KMATRIX_ASL]"
		"Pardiso done. \n");

  /* v3 */
  /* */
  somefile = fopen("result.txt", "w");
  for(i=0; i<K_nrows; i++){
    fprintf(somefile, "\t%d", i);
		for(j=0; j<n_dof; j++){
    	fprintf(somefile, "\t%f", *(x_+ j * K_nrows + i));
    }
      fprintf(somefile, "\n");
      }
  fclose(somefile);
  somefile = fopen("result_unscaled.txt", "w");
  for(i=0; i<K_nrows; i++){
    /*fprintf(somefile, "\t%d", i);*/
		for(j=0; j<n_dof; j++){
			*(x_+ j * K_nrows + i) = *(x_+ j * K_nrows + i) * exp(S_scale[i]);
    	fprintf(somefile, "\t%f", *(x_+ j * K_nrows + i));
    }
      fprintf(somefile, "\n");
      }
  fclose(somefile);

  memset(positions_rh, 0, sizeof(int)*n_var);

  for(i=0; i<n_dof; i++){
  	j = hr_point[i];
  	positions_rh[j] = i;
  }
  somefile = fopen("result_red_hess.txt", "w");
  for(i=0; i<n_dof; i++){
		for(j=0; j<n_dof; j++){
    	fprintf(somefile, "\t%f", *(x_+ j * K_nrows + hr_point[i]));
    }
    fprintf(somefile, "\n");
  }
  fclose(somefile);


  somefile = fopen("result_primal_dual.txt", "w");
  
  for(i=0; i<K_nrows; i++){
    fprintf(somefile, "\t%f\n", s_star[i]);
  }
  fclose(somefile);

  suf_iput(reg_suffix_name[1], ASL_Sufkind_var, positions_rh);


  write_sol(ter_msg, s_star, s_star + n_var, 0);
  free(positions_rh);
	free(reg_suffix_name[0]);
	free(reg_suffix_name[1]);
	free(reg_suffix_name);
	free(nz_row_w);
	free(nz_row_a);
	free(md_off_w);
	free(rhs_baksolve);
	free(x_);
	free(hr_point);
	free(dp_);
	free(s_star);
	/* suf_name = (char **)malloc(sizeof(char *)*n_rhs); */
	free(S_scale);
	free(c_flag);

	

	for(i=0; i<n_rhs; i++){
		free(rhs_name[i]);
	}
	free(rhs_name);

	ASL_free(&asl);
	free(suf_ptr);
	free(Arow);
	free(Acol);
	free(Aij);
	free(Wrow);
	free(Wcol);
	free(Wij);
	free(Wr_t);
	free(Wc_t);
	free(Wi_t);
	free(Krow);
	free(Kcol);
	free(Kij);
	free(Kr_strt);
	free(rhs_ptr);
	return 0;
}
