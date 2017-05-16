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
#include "sort_by_column.h"
#include "w_append_nz.h"
#include "k_assemble_cc.h"
#include "kmalloc.h"
#include "mc30_driver.h"
#include "pardiso_driver.h"

#include "get_jac_asl.h"
#include "get_hess_asl.h"

#include "find_inequalities.h"
#include "assemble_corrector_rhs.h"
/*#include "assemble_rhsds.h" */
#include "assemble_rhsdsv3.h"

static int dumm = 1;
static I_Known dumm_kw = {2, &dumm};
static int n_rhs = 0;
static int l_over = 0;
static I_Known l_over_kw = {1, &l_over};

/* keywords */
static keyword keywds[] = {
  KW("smth", IK_val, &dumm_kw, "Cheers mate the cavalry is here"),
  KW("n_rhs", I_val, &n_rhs, "Number or right hand sides"),	
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
	static char sword[] = "rhs_";
	SufDesc *var_f=NULL;
	SufDesc *con_f=NULL;
	SufDesc **rhs_ptr=NULL;
	SufDecl *tptr=NULL;

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
	int *c_flag=NULL;
	char **suf_name=NULL;

	/* the pointer to array that will contain */
	/* the rhs */
	real *rhs_baksolve = NULL;

  FILE *somefile;
  fint nerror;
  fint nzA;

  char ter_msg[] = {"I[KMATRIX]...\t[KMATRIX_ASL]"
	"All done it seems."};

	/* The memory allocation */
	asl = ASL_alloc(ASL_read_pfgh);


	s = getstops(argv, &Oinfo);
	/*s = getstub(&argv,&Oinfo); */


	if (!s) {
		return 1;
	}
	else {
		printf("I[KMATRIX]...\t[KMATRIX_ASL]"
			"File read succesfull\n");
				}

	if (n_rhs == 0){
		printf("E[KMATRIX]...\t[KMATRIX_ASL]"
			"No n_rhs declared\n");
		return -1;
	}

	printf("I[KMATRIX]...\t[KMATRIX_ASL]"
	"Number of Right hand sides %d\n", n_rhs);
/*	printf("Number of Right hand sides %d\n", n_rhs); */

	/* Declare suffix array pointer */
	tptr = (SufDecl *)malloc(sizeof(SufDecl)*(n_rhs+2));

	/* First two suffixes */
	tptr->name = "var_flag";
	tptr->table = 0;
	tptr->kind = ASL_Sufkind_var;
	tptr->nextra = 0;

	(tptr + 1)->name = "con_flag";
	(tptr + 1)->table = 0;
	(tptr + 1)->kind = ASL_Sufkind_con;
	(tptr + 1)->nextra = 0;

	/* Suffix for the right hand side creation */
	suf_name = (char **)malloc(sizeof(char *)*n_rhs);

	/* this bit writes names of the rhs suffixes */
	for(i=0; i < n_rhs; i++){
		char numeric_rhs[5];
		numeric_rhs[0] = '\0';
	  sprintf(numeric_rhs, "%d", i);
	  suf_name[i] = (char *)malloc(sizeof(char)*10);
	  *(suf_name[i]) = '\0';
	  strcat(suf_name[i], sword);
	  strcat(suf_name[i], numeric_rhs);
/*	  printf("The name of the suffix %s\n", suf_name[i]); */
	}

	for(i=0; i<n_rhs; i++){
		printf("NAME %s\n", suf_name[i]);
	}

	for(i=0; i<n_rhs; i++){
		(tptr + i + 2)->name = suf_name[i];
		(tptr + i + 2)->table = 0;
		(tptr + i + 2)->kind = ASL_Sufkind_con|ASL_Sufkind_real;
		(tptr + i + 2)->nextra = 0;
	}

	/* Declare suffixes */
	suf_declare(tptr, (n_rhs + 2));

	/* dhis bit setups ASL components e.g. n_var, n_con, etc. */
	f = jac0dim(s, (fint)strlen(s));

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
			free(suf_name[i]);
		}
		free(suf_name);
		ASL_free(&asl);
		free(tptr);
		return -1;
	}


	con_f = suf_get("con_flag", ASL_Sufkind_con);
	var_f = suf_get("var_flag", ASL_Sufkind_var);

	if(var_f->u.i == NULL){
    fprintf(stderr, "E[KMATRIX]...\t[KMATRIX_ASL]"
    	"u.i empty, no var_flag declared.\n");
    /*return -1; */
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
   *(rhs_ptr + i)= suf_get(suf_name[i], ASL_Sufkind_con|ASL_Sufkind_real);
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
	get_jac_asl(asl,x,Acol,Arow,Aij,nzc,&nerror);
	get_hess_asl(asl,x,&Wcol,&Wrow,&Wij, n_var, n_con, n_obj, &nnzw, lambda, &nerror);
	/* */
	printf("I[KMATRIX]...\t[KMATRIX_ASL]"
	"Nonzeroes in the sparse hessian %d\n", nnzw);
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
	sortcol(Arow, Acol, Aij, nzA, n_var);
	f_debug = fopen("somefile.txt", "w");
	for(i=0; i < nzA; i++){
    fprintf(f_debug, "%d\t%d\t%.g\n", Arow[i], Acol[i], Aij[i]);
	}
	fclose(f_debug);


	/*wnzappnd(Wrow, Wcol, Wij, nnzw, n_var); */
	/*w_col_sort(Wrow, Wcol, Wij, nnzw, n_var); */
	/* calculate the space required for the KKT_matrix */
	k_space = k_malloc(Wrow, Wcol, Wij, nnzw, n_var, n_con, nzA);

	/* temporal space for w_appended */
	wa_tsp = k_space - n_con - nzA;
	Krow = (fint *)malloc(sizeof(fint) * k_space);
	Kcol = (fint *)malloc(sizeof(fint) * k_space);
	Kij  = (real *)malloc(sizeof(real) * k_space);
	/* wnzappnd(fint *Wrow, fint *Wcol, real *Wij, fint Wnz, fint nvar, */
	/* fint *Wr_new, fint *Wc_new, real *Wi_new, fint Wnz_new) */

	Wr_t = (fint *)malloc(sizeof(fint) * wa_tsp);
	Wc_t = (fint *)malloc(sizeof(fint) * wa_tsp);
	Wi_t = (real *)malloc(sizeof(real) * wa_tsp);

	/* append 0's to the main diagonal as necessary */
	wnzappnd(Wrow, Wcol, Wij, nnzw, n_var, Wr_t, Wc_t, Wi_t, wa_tsp);
	/* fint *Wrow, fint *Wcol, real *Wij, fint Wnz, fint nvar, fint ncon, */
	/* fint *Arow, fint *Acol, real *Aij, fint Anz */
	f_debug = fopen("appended_hessian.txt", "w");
	for(i=0; i<wa_tsp; i++){
    fprintf(f_debug, "\t%ld\t%ld\t%.g\n", Wr_t[i], Wc_t[i], Wi_t[i]);
	}
	fclose(f_debug);


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
  rhs_baksolve = (real *)calloc(K_nrows * n_rhs, sizeof(real));
  x_           = (real *)calloc(K_nrows * n_rhs, sizeof(real));
  dp_					 = (real *)calloc(n_rhs, sizeof(real)); /* array that contains delta_p */

  
	/* */
	assemble_rhsds(n_rhs, K_nrows, rhs_baksolve, dp_, n_var, n_con, rhs_ptr);
  
	for(i=0; i<n_rhs; i++){
		printf("rhs_%d, delta_p %f\n", i, dp_[i] );
	}

  /* scale matrix */
  for(i=0; i< k_space; i++){
  	Kij[i] = Kij[i] * exp(S_scale[Kcol[i]-1]) * exp(S_scale[Krow[i]-1]);
  }
  
  /* scale rhs by [n][rhs] */
  /*
  for(i=0; i< K_nrows; i++){
  	for(j=0; j < n_rhs; j++)
  		*(*(rhs_baksolve+i)+j) = *(*(rhs_baksolve + i)+j) * exp(S_scale[i]);
  }
  */

 	/*assert(&(rhs_baksolve[300][1]) == (*(rhs_baksolve + 300)+1)); */
  /*
  // scale rhs by [rhs][n] 
  for(i=0; i< n_rhs; i++){
  	for(j=0; j < K_nrows; j++)
  		*(*(rhs_baksolve+i)+j) = *(*(rhs_baksolve + i)+j) * exp(S_scale[j]);
  }
	*/

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
	pardiso_driver(Kr_strt, Kcol, Kij, K_nrows, k_space, n_rhs, rhs_baksolve, x_); 
      
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
    fprintf(somefile, "\t%d", i);
		for(j=0; j<n_rhs; j++){
			*(x_+ j * K_nrows + i) = *(x_+ j * K_nrows + i) / exp(S_scale[i]);
    	fprintf(somefile, "\t%f", *(x_+ j * K_nrows + i));
    }
      fprintf(somefile, "\n");
      }
  fclose(somefile);
  write_sol(ter_msg, x, lambda, 0);
/* */

  /*
  somefile = fopen("result.txt", "w");
  for(i=0; i<K_nrows; i++){
    fprintf(somefile, "\t%d", i);
		for(j=0; j<n_rhs; j++){
    	fprintf(somefile, "\t%f", *(*(x_+j)+i));
    }
      fprintf(somefile, "\n");
      }
  fclose(somefile);
*/
/*
  somefile = fopen("result.txt", "w");
  for(i=0; i<K_nrows; i++){
    fprintf(somefile, "\t%d", i);
   for(j=0; j<n_rhs; j++){
    fprintf(somefile, "\t%f", *(*(x_+i)+j));
    }
  fprintf(somefile, "\n");
  }
  fclose(somefile);

  
	for(i=0; i<n_rhs; i++){ // by [rhs][n]
	//for(i=0; i<K_nrows; i++){ // by [n][rhs]
		free(*(rhs_baksolve+i));
		free(*(x_+i));
	}
  */
	free(rhs_baksolve);
	free(x_);
	free(dp_);
	/* suf_name = (char **)malloc(sizeof(char *)*n_rhs); */
	free(S_scale);
	free(c_flag);

	

	for(i=0; i<n_rhs; i++){
		free(suf_name[i]);
	}
	free(suf_name);

	ASL_free(&asl);
	free(tptr);
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
