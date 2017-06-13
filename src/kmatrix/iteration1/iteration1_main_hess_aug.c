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
/*#include "sort_by_column.h" */
#include "w_append_nz.h"
#include "k_assemble_cc.h"
#include "kmalloc.h"
#include "mc30_driver.h"
#include "pardiso_driver.h"

#include "get_jac_asl_aug.h"
#include "get_hess_asl_aug.h"

#include "find_inequalities.h"
#include "assemble_corrector_rhs.h"
/*#include "assemble_rhsds.h" */
#include "assemble_rhsd_red_hess.h"
#include "sens_update_driver.h"
#include "suffix_decl_hand.h"

static int dumm = 1;
static I_Known dumm_kw = {2, &dumm};
static int n_rhs = 0;
static int n_dof = 0;
static int l_over = 0;
static I_Known l_over_kw = {1, &l_over};

/* keywords */
static keyword keywds[] = {
  KW("smth", IK_val, &dumm_kw, "Cheers mate the cavalry is here"),
  KW("n_dof", I_val, &n_dof, "another keyword"),
  KW("n_rhs", I_val, &n_rhs, "Number of right hand sides"),	
  KW("no_lambda", IK_val, &l_over_kw, "Override multiplier check(lambda)")
};
static char banner[] = {"[KMATRIX] written by DT\n\n"};
static Option_Info Oinfo = 
{"k_matrix", banner, "k_matrix_options", keywds, nkeywds};
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

	fint *Acol=NULL, *Arow=NULL;
	real *Aij=NULL;
	fint *Wcol=NULL, *Wrow=NULL;
	real *Wij=NULL;
	fint *Kcol=NULL, *Krow=NULL;
	real *Kij=NULL;
	fint *Kr_strt=NULL;
	fint K_nrows;

	real *S_scale=NULL;
	fint k_space;

	fint *Wc_t=NULL, *Wr_t=NULL;
	real *Wi_t=NULL;
	fint wa_tsp;

	real *x_=NULL;
	real *dp_=NULL;
	/*real *gf= NULL; */
	real *s_star = NULL;
	int *c_flag=NULL;
	char **rhs_name=NULL;
	int *nz_row_w=NULL;
	int *nz_row_a=NULL;
	int *md_off_w=NULL;
	int miss_nz_w;

	/* the pointer to array that will contain */
	/* the rhs */
	real *rhs_baksolve = NULL;

  FILE *somefile;
  fint nerror;
  fint nzA;
  int *hr_point = NULL;

  char ter_msg[] = {"I[KMATRIX]...\t[KMATRIX_ASL]"
	"All done it seems."};

	printf("I[PRE-TEST]...\t[KMATRIX_ASL]"
	"Number of Right hand sides %d\n", n_rhs);

	/*if(n_rhs > 0){
		suf_ptr = (SufDecl *)malloc(sizeof(SufDecl)*(n_rhs+2));
		rhs_name = (char **)malloc(sizeof(char *)*n_rhs);
		for(i=0; i<n_rhs; i++){
			rhs_name[i] = (char *)malloc(sizeof(char) * 32); 
		}
		suffix_decl_hand(n_rhs, suf_ptr, rhs_name);
	}
	else{
		suf_ptr = (SufDecl *)malloc(sizeof(SufDecl)*2);
		suffix_decl_hand(n_rhs, suf_ptr, rhs_name);
	}
	*/


	/* The memory allocation */
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
		printf("E[KMATRIX]...\t[KMATRIX_ASL]"
			"No n_rhs declared\n");
		/*exit(-1);*/
	}


	if(n_rhs > 0){
		suf_ptr = (SufDecl *)malloc(sizeof(SufDecl)*(n_rhs+2));
		rhs_name = (char **)malloc(sizeof(char *)*n_rhs);
		for(i=0; i<n_rhs; i++){
			rhs_name[i] = (char *)malloc(sizeof(char) * 32); /* 32 bit long digit; why not?*/
		}
		suffix_decl_hand(n_rhs, suf_ptr, rhs_name);
	}
	else{
		suf_ptr = (SufDecl *)malloc(sizeof(SufDecl)*2);
		suffix_decl_hand(n_rhs, suf_ptr, rhs_name);
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

	if (n_dof > (n_var - n_con) || n_dof == 0){
		printf("E[KMATRIX]...\t[KMATRIX_ASL]"
			"No invalid number of n_dof declared\n");
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


	con_f = suf_get("con_flag", ASL_Sufkind_con);
	var_f = suf_get("dof_v", ASL_Sufkind_var);

	if(var_f->u.r == NULL && var_f->u.i == NULL){
    fprintf(stderr, "E[KMATRIX]...\t[KMATRIX_ASL]"
    	"u.i empty, no var_flag declared.\n");
    return -1;
	}
	if(con_f->u.i == NULL){
		fprintf(stderr, "E[KMATRIX]...\t[KMATRIX_ASL]"
			"u.i empty, no con_flag declared.\n");
	/*	return -1; */
	}

	/*SufDesc *smptr0, *smptr1; */

	/* smptr0 = suf_get("rhs_0", ASL_Sufkind_con); */
	/* smptr1 = suf_get("rhs_1", ASL_Sufkind_con); */
	rhs_ptr = (SufDesc **)malloc(sizeof(SufDesc *) * n_rhs);

	for(i=0; i < n_rhs; i++){
   *(rhs_ptr + i)= suf_get(rhs_name[i], ASL_Sufkind_con|ASL_Sufkind_real);
  	if((*(rhs_ptr + i))->u.r == NULL){
		  fprintf(stderr, "E[KMATRIX]...\t[KMATRIX_ASL]"
		  	"u.r empty, no rhs values declared for rhs_%d.\n", i);
		  return -1;
  	}
	}
	
	/*
	real temp;	
	for(i=0; i < n_con; i++){
		printf("%d", i);
		for(int k=0; k < n_rhs; k++){
			temp = (*(rhs_ptr + k))->u.r[i];
			printf("\t%f", temp);
		}
		printf("\n");
	}
	*/

	c_flag = (int *)malloc(sizeof(int) * n_con);
	
	/*constraints flags */
	find_ineq_con(n_con, LUrhs, c_flag);

	/* Row and colum for the triplet format A matrix */
	/* size of the number of nonzeroes in the constraint jacobian */
	Acol = (fint *)malloc(sizeof(fint)*nzc);
	Arow = (fint *)malloc(sizeof(fint)*nzc);
	Aij  = (real *)malloc(sizeof(real)*nzc);

	
	nerror = 0;
	
	/* get_jac_asl function */
	/* is it better to pass asl as a reference? */
	/*get_grad_f(asl, x, n_var, n_obj, gf, &nerror); */
	get_jac_asl_aug (asl, x, Acol, Arow, Aij, n_var, n_con, nzc, &nerror, &nz_row_a);
	get_hess_asl_aug(asl, x, &Wcol, &Wrow, &Wij, n_var, n_con, n_obj,
	 &nnzw, lambda, &nerror, &nz_row_w, &md_off_w, miss_nz_w);
	/* */
	exit(-1);
	/*printf("nonzeroes in the sparse hessian %d\n", nnzw); */
	/*printf("print dummy %d\n", dumm); */
	
	if (dumm == 2){
    printf("I[KMATRIX]...\t[KMATRIX_ASL]"
    "Cheers mate, the cavalry is here!");
	}

	/*printf("size of s %d\n", sizeof(f)); */
	printf("I[KMATRIX]...\t[KMATRIX_ASL]"
		"Number of Right hand sides %d\n", n_rhs);
	printf("I[KMATRIX]...\t[KMATRIX_ASL]"
		"Number of variables %d\n", n_var);
	printf("I[KMATRIX]...\t[KMATRIX_ASL]"
		"Number of constraints %d\n", n_con);
	
	nzA = nzc;
	/* Reorder the jacobian */
	/* Done by now */
	


	/*wnzappnd(Wrow, Wcol, Wij, nnzw, n_var); */
	/*w_col_sort(Wrow, Wcol, Wij, nnzw, n_var); */
	/* calculate the space required for the KKT_matrix */
	
	k_space = nzA + nnzw + miss_nz_w + n_var;

	/* temporal space for w_appended */
	wa_tsp = k_space - n_con - nzA;
	Krow = (fint *)malloc(sizeof(fint) * k_space);
	Kcol = (fint *)malloc(sizeof(fint) * k_space);
	Kij  = (real *)malloc(sizeof(real) * k_space);
	

	/* the question is if we do not append nz to the W matrix, can k_assemble  */
	/* still do its work? */
	/* appearently, yes it will */

	/* row starts */
	K_nrows = n_var + n_con;
	/* number of rows of the KKT matrix */
	printf("I[KMATRIX]...\t[KMATRIX_ASL]"
		"n_rows KKT %d\n", K_nrows);

	Kr_strt = (fint *)malloc(sizeof(fint) * (K_nrows + 1));
	S_scale = (real *)malloc(sizeof(real) * K_nrows);

	
	/*real * Crhs=NULL; */
	/* assemble_corrector_rhs(asl, x, lambda, n_var, n_con, Arow, Acol, Aij, nzA, Crhs, &nerror, c_flag); */

	k_assemble_cc(Wr_t, Wc_t, Wi_t, wa_tsp, n_var, n_con, 
		Arow, Acol, Aij, nzA, 
		Krow, Kcol, Kij, k_space,
	        Kr_strt);

	mc30driver(K_nrows, k_space, Kij, Krow, Kcol, S_scale);

	/* by [n][rhs] */
  /*
  rhs_baksolve = (real **)malloc(sizeof(real *) * K_nrows);
	x_ =           (real **)malloc(sizeof(real *) * K_nrows);
	assert(rhs_baksolve != NULL);
	assert(x_ != NULL);
	for(i = 0; i < K_nrows; i++){
		*(rhs_baksolve + i) = (real *)calloc(n_rhs, sizeof(real));
		*(x_ + i) =           (real *)calloc(n_rhs, sizeof(real));
		assert(*(rhs_baksolve + i) != NULL);
		assert(*(x_ + i) != NULL);
	}
  */

	/* */
  /* by [rhs][n] */
  rhs_baksolve = (real *)calloc(K_nrows * (n_dof), sizeof(real));
  x_           = (real *)calloc(K_nrows * (n_dof), sizeof(real));
  hr_point		 = (int *)malloc(n_dof * sizeof(int));

  s_star			 = (real *)calloc(K_nrows, sizeof(real)); /* array that contains the primal dual update */

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
  printf("Here0!\n");
  assemble_rhs_red_hess(rhs_baksolve, n_var, n_con, n_dof, var_f, hr_point);
  printf("Here1!\n");
  
	
  /* scale matrix */
  for(i=0; i< k_space; i++){
  	Kij[i] = Kij[i] * exp(S_scale[Kcol[i]-1]) * exp(S_scale[Krow[i]-1]);
  }
  
  
  /* v3 */
  for(i=0; i< n_rhs; i++){
  	for(j=0; j < K_nrows; j++){
  		*(rhs_baksolve + i*K_nrows + j) = *(rhs_baksolve + i*K_nrows + j) * exp(S_scale[j]);
  	}
  }

	somefile = fopen("rhs_sens_scaled", "w");

	/*
 	for(j=0; j < K_nrows; j++){
 		for(i=0; i < n_rhs; i++){
 			fprintf(somefile, "\t%f\t", *(*(rhs_baksolve + i)+j) );
 		}
 		fprintf(somefile, "\n");
 	}
 	fclose(somefile); */

 	/* v3 */
 	for(j=0; j < K_nrows; j++){
 		for(i=0; i < n_rhs; i++){
 			fprintf(somefile, "\t%f\t", *(rhs_baksolve + i*K_nrows + j) );
 		}
 		fprintf(somefile, "\n");
 	}
 	fclose(somefile);

 
  /*  factorize the matrix */
	pardiso_driver(Kr_strt, Kcol, Kij, K_nrows, k_space, n_dof, rhs_baksolve, x_);
      
  printf("I[KMATRIX]...\t[KMATRIX_ASL]"
		"Pardiso done. \n");

  /* v3 */
  /* */
  somefile = fopen("result.txt", "w");
  for(i=0; i<K_nrows; i++){
    fprintf(somefile, "\t%d", i);
		for(j=0; j<n_rhs; j++){
    	fprintf(somefile, "\t%f", *(x_+ j * K_nrows + i));
    }
      fprintf(somefile, "\n");
      }
  fclose(somefile);
  somefile = fopen("result_unscaled.txt", "w");
  for(i=0; i<K_nrows; i++){
    /*fprintf(somefile, "\t%d", i);*/
		for(j=0; j<n_dof; j++){
			*(x_+ j * K_nrows + i) = *(x_+ j * K_nrows + i) / exp(S_scale[i]);
    	fprintf(somefile, "\t%f", *(x_+ j * K_nrows + i));
    }
      fprintf(somefile, "\n");
      }
  fclose(somefile);


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


  write_sol(ter_msg, s_star, s_star + n_var, 0);
  
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
