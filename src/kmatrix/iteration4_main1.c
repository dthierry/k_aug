/* @source kmatrix_b0.c
** beta 0
** April 18th, 2017
** @author: David Thierry (dmolinat@andrew.cmu) dav0@lb2016-1

********************************************************************************

@main ********************************************
**
** Reads nl file, allocates data structures, calls assembling funcs
** ToDo:
** Need to implement rhs and red hess in a single program
** Write program that takes suffixes from dot_prod calculation and performs the
** Sensitivity step. 
** @param [r] stub
** @param [r] KWs
** @@
*******************************************************************************/
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "getstub.h"
#include "mc30_driver.h"
#include "pardiso_driver.h"
#include "get_jac_asl_aug.h"
#include "get_hess_asl_aug.h"
#include "find_inequalities.h"
#include "assemble_rhs_rh.h"
#include "suffix_decl_hand.h"
#include "csr_driver.h"
#include "sigma_compute.h"
#include "mu_adjust_primal.h"
#include "dsyev_driver.h"

#define NUM_REG_SUF 4

static real not_zero = 1.84e-04;
static int dumm = 1;
static I_Known dumm_kw = {2, &dumm};
static int n_rhs = 0;
static int l_over = 0;
static I_Known l_over_kw = {1, &l_over};

static char _dot_pr_f[] = {"dot_prod"};
static char name1[] = {"smth"};
static char _e_eval[] = {"eig_rh"};
static char _n_rhsopt_[] = {"n_rhs"};
static char _no_barrieropt_[] = {"no_barrier"};
static char _no_lambdaopt_[] = {"no_lambda"};
static char _no_scaleopt_[]  = {"no_scale"};
static char _not_zero[] = {"not_zero"};



static int dot_prod_f = 0;
static I_Known dot_p_kw = {1, &dot_prod_f};

static int eig_rh_eval = 0;
static I_Known e_eval_kw = {1, &eig_rh_eval};

static int no_barrier = 1;
static I_Known nbarrier_kw = {0, &no_barrier};

static int no_scale = 1;
static I_Known nscale_kw = {0, &no_scale};


/*static char dof_v[] = {"dof_v"};*/

/* keywords, they must be in alphabetical order! */
static keyword keywds[] = {
	KW(_dot_pr_f, IK_val, &dot_p_kw, _dot_pr_f),
	KW(_e_eval, IK_val, &e_eval_kw, _e_eval),
  KW(_n_rhsopt_ , I_val, &n_rhs, _n_rhsopt_),	
  KW(_no_barrieropt_ , IK_val, &nbarrier_kw, _no_barrieropt_),
  KW(_no_lambdaopt_ , IK_val, &l_over_kw, _no_lambdaopt_),  
  KW(_no_scaleopt_ , IK_val, &nscale_kw, _no_scaleopt_),  
  KW(_not_zero , D_val, &not_zero, _not_zero),  
  KW(name1 , IK_val, &dumm_kw, name1),
};
static char banner[] = {"[KMATRIX] written by DT\n\n"};
static char _k_[] = {"K_augmented"};
static char _k_o_[] = {"K_augmented_options"};
static Option_Info Oinfo;




int main(int argc, char **argv){
	ASL *asl;
	FILE *f;

	/* SufDesc *some_suffix; */
	int i, j;
	int n_dof=0;
	int nnzw; /* let's try this */
	real *x=NULL, *lambda=NULL;
	char *s=NULL;
	SufDesc *var_f=NULL;

	SufDesc *suf_zL = NULL;
	SufDesc	*suf_zU = NULL;
	real *z_L=NULL, *z_U=NULL, *sigma=NULL;

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
	fint nzK;

	fint *Wc_t=NULL, *Wr_t=NULL;
	real *Wi_t=NULL;


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

  char ter_msg[] = {"I[KMATRIX]...[KMATRIX_ASL]"
	"All done it seems."};

	unsigned n_r_suff = NUM_REG_SUF;
	/* Suffix names; yes, I know. */
	char _sfx_1[] = {"dof_v"};
	char _sfx_2[] = {"rh_name"};
	char _sfx_3[] = {"ipopt_zL_in"};
	char _sfx_4[] = {"ipopt_zU_in"};
	


	Oinfo.sname = _k_;
	Oinfo.bsname = banner;
	Oinfo.opname = _k_o_;
	Oinfo.keywds = keywds;
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
	
	if (l_over){
		printf("W[KMATRIX]...\t[KMATRIX_ASL]"
			"Multiplier check override.\n");
	}
	if (dot_prod_f){
		printf("W[KMATRIX]...\t[KMATRIX_ASL]"
			"Dot product preparation.\n");
	}


	/* Allocate suffix names (regular)*/
	reg_suffix_name =   (char **)malloc(sizeof(char *) * n_r_suff);
	for(i=0; i < (int)n_r_suff; i++){
		reg_suffix_name[i] = (char *)malloc(sizeof(char) * 16 );
		reg_suffix_name[i][0] = '\0';
	}
	
	strcat(reg_suffix_name[0], _sfx_1);
	strcat(reg_suffix_name[1], _sfx_2);
	strcat(reg_suffix_name[2], _sfx_3);
	strcat(reg_suffix_name[3], _sfx_4);

	if(n_rhs > 0){
		suf_ptr = (SufDecl *)malloc(sizeof(SufDecl)*(n_rhs + n_r_suff));
		rhs_name = (char **)malloc(sizeof(char *)*n_rhs);
		for(i=0; i<n_rhs; i++){
			rhs_name[i] = (char *)malloc(sizeof(char) * 32); /* 32 bit long digit;
			 why not?*/
		}
		suffix_decl_hand(suf_ptr, reg_suffix_name, rhs_name, n_r_suff, n_rhs);
	}
	else{
		suf_ptr = (SufDecl *)malloc(sizeof(SufDecl) * n_r_suff);
		suffix_decl_hand(suf_ptr, reg_suffix_name, rhs_name, n_r_suff, n_rhs);
	}


	printf("I[KMATRIX]...\t[KMATRIX_ASL]"
	"Number of Right hand sides %d\n", n_rhs);
	
	

	/* Declare suffixes */
	if(n_rhs > 0){suf_declare(suf_ptr, (n_rhs + n_r_suff));}
	else{suf_declare(suf_ptr, n_r_suff);}
	

	/* dhis bit setups ASL components e.g. n_var, n_con, etc. */
	f = jac0dim(s, (fint)strlen(s));


	printf("I[KMATRIX]...\t[KMATRIX_ASL]"
		"Number of Right hand sides: %d\n", n_rhs);
	printf("I[KMATRIX]...\t[KMATRIX_ASL]"
		"Number of variables       : %d\n", n_var);
	printf("I[KMATRIX]...\t[KMATRIX_ASL]"
		"Number of constraints     : %d\n", n_con);
	printf("I[KMATRIX]...\t[KMATRIX_ASL]"
		"Number of valid n_dof     : %d\n", n_var - n_con );


	if ((n_var - n_con) < 0){
		printf("E[KMATRIX]...\t[KMATRIX_ASL]"
			"nvar < ncon. This problem is not valid.\n");
		exit(-1);
	}
	
	x 		 = X0  = M1alloc(n_var*sizeof(real));
	lambda = pi0 = M1alloc(n_con*sizeof(real));
	obj_no = 0;
	/* need to do part of changing sign for y */

	pfgh_read(f, ASL_findgroups);

	
	/* NEED TO FIX THIS	*/
	if(lambda==NULL && l_over == 0){
		printf("E[KMATRIX]...\t[KMATRIX_ASL]"
	"Constraint Multipliers not declared(suffix dual), abort\n");
		for(i=0; i < (int)n_r_suff; i++){free(reg_suffix_name[i]);}
		free(reg_suffix_name);
		if(n_rhs){
			for(i=0; i<n_rhs; i++){
				free(rhs_name[i]);
			}
			free(rhs_name);
		}
		free(suf_ptr);
		ASL_free(&asl);
		exit(-1);
	}

	var_f = suf_get(reg_suffix_name[0], ASL_Sufkind_var);
	
	suf_zL = suf_get(reg_suffix_name[2], ASL_Sufkind_var| ASL_Sufkind_real); 
	suf_zU = suf_get(reg_suffix_name[3], ASL_Sufkind_var| ASL_Sufkind_real); 

	z_L = (real *)malloc(sizeof(real) * n_var);
	z_U = (real *)malloc(sizeof(real) * n_var);

	memset(z_L, 0, sizeof(real) * n_var);
	memset(z_U, 0, sizeof(real) * n_var);


	if(!(suf_zL->u.r)){
		fprintf(stderr, "W[KMATRIX_ASL]...\t[KMATRIX_ASL]"
    	"No ipopt_zL_out suffix declared, setting zL = 0.\n");
	}
	else{
		for(i=0; i< n_var; i++){
			z_L[i] = suf_zL->u.r[i];
		}
	}
	if(!(suf_zU->u.r)){
		fprintf(stderr, "W[KMATRIX_ASL]...\t[KMATRIX_ASL]"
    	"No ipopt_zU_out suffix declared, setting zU = 0.\n");
	}
	else{
		for(i=0; i< n_var; i++){
			z_U[i] = suf_zU->u.r[i];
		}
	}


	somefile = fopen("primal0.txt", "w");
	for(i=0; i< n_var; i++){
		fprintf(somefile, "%.g\n", x[i]);
	}
	fclose(somefile);

	mu_adjust_x(n_var, x, LUv, z_L, z_U);
	
	somefile = fopen("primal1.txt", "w");
	for(i=0; i< n_var; i++){
		fprintf(somefile, "%.g\n", x[i]);
	}
	fclose(somefile);

	sigma = (real *)malloc(sizeof(real) * n_var);
	memset(sigma, 0, sizeof(real) * n_var);

	if(var_f->u.r == NULL && var_f->u.i == NULL){
    fprintf(stderr, "E[KMATRIX]...\t[KMATRIX_ASL]"
    	"suffix empty no n_dof declared!\n");
    exit(-1);
	}

	compute_sigma(asl, n_var, x, suf_zL, suf_zU, z_L, z_U, sigma, not_zero);
	
	
	/* Is this gonna work? */
	if(n_rhs > 1){
		rhs_ptr = (SufDesc **)malloc(sizeof(SufDesc *) * n_rhs);
		for(i=0; i < n_rhs; i++){
	   *(rhs_ptr + i)= suf_get(rhs_name[i], ASL_Sufkind_con|ASL_Sufkind_real);
	  	if((*(rhs_ptr + i))->u.r == NULL){
			  fprintf(stderr, "E[KMATRIX]...\t[KMATRIX_ASL]"
			  	"No rhs values declared for rhs_%d.\n", i);
			  exit(-1);
	  	}
		}
	}
	
	c_flag = (int *)malloc(sizeof(int) * n_con); /* Flags for ineq or equalities*/
	
	/*constraints flags */
	find_ineq_con(n_con, LUrhs, c_flag); /* Find the inequality constraints */

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

	if(no_barrier){
	/* Add barrier term to the main diagonal */
		for(i=0; i<n_var; i++){
			j = md_off_w[i];
			Wij[j] += sigma[i];
		}
		printf("I[KMATRIX]...\t[KMATRIX_ASL]"
		"Barrier term added.\n");
	}
	

	nzA = nzc; 
	nzW = nnzw + miss_nz_w; /* Recomputed number of nz in the Hessian-of-Lagrangian */


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
  assemble_rhs_rh(&rhs_baksolve, n_var, n_con, &n_dof, var_f, &hr_point);

  x_           = (real *)calloc(K_nrows * (n_dof), sizeof(real));
  positions_rh = (int *)malloc(n_var * sizeof(int));
  /*for(i=0; i<n_dof; i++){
  	printf("i %d, hr %d\n", i, hr_point[i]);
  }*/
  
	
  /* scale matrix & rhs*/
  if(no_scale > 0){
  	for(i=0; i< nzK; i++){
  		Kij[i] = Kij[i] * exp(S_scale[Kcol[i]-1]) * exp(S_scale[Krow[i]-1]);
  	}
  	for(i=0; i< n_dof; i++){
	  	for(j=0; j < K_nrows; j++){
	  		*(rhs_baksolve + i*K_nrows + j) = 
	  		*(rhs_baksolve + i*K_nrows + j) * exp(S_scale[j]);
	  	}
  	}
  	somefile = fopen("rhs_sens_scaled", "w");
 		for(j=0; j < K_nrows; j++){
 			for(i=0; i < n_dof; i++){
 				fprintf(somefile, "\t%f\t", *(rhs_baksolve + i*K_nrows + j) );
 			}
 			fprintf(somefile, "\n");
 		}
 		fclose(somefile);
  }
  else{
  	printf("W[KMATRIX]...\t[KMATRIX_ASL]"
			"The scaling has been skipped. \n");
	}
 
  /* factorize the matrix */
	pardiso_driver(Kr_strt, Kcol, Kij, K_nrows, nzK, n_dof, rhs_baksolve, x_);
      
  printf("I[KMATRIX]...\t[KMATRIX_ASL]"
		"Pardiso done. \n");

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

  /*somefile = fopen("result_red_hess_scaled.txt", "w");
  for(i=0; i<n_dof; i++){
		for(j=0; j<n_dof; j++){
			fprintf(somefile, "\t%.g", (*(x_+ j * K_nrows + hr_point[i]) + *(x_+ i * K_nrows + hr_point[j])) * 0.5 );
    }
    fprintf(somefile, "\n");
  }
  fclose(somefile);*/


  somefile = fopen("result_unscaled.txt", "w");
  if(no_scale > 0){
  	for(i=0; i<K_nrows; i++){
			for(j=0; j<n_dof; j++){
				*(x_+ j * K_nrows + i) = *(x_+ j * K_nrows + i) * exp(S_scale[i]);
    		fprintf(somefile, "\t%f", *(x_+ j * K_nrows + i));
    	}
      fprintf(somefile, "\n");
    }  	
  }
  fclose(somefile);

  somefile = fopen("dot_in.in", "w");
  for(i=0; i<n_dof; i++){
		for(j=0; j<K_nrows; j++){
    	fprintf(somefile, "\t%.g\n", *(x_+ i * K_nrows + j));
    }
      }
  fclose(somefile);


  memset(positions_rh, 0, sizeof(int)*n_var);
  somefile = fopen("sigma_warnings.txt", "w");
  for(i=0; i<n_dof; i++){
  	j = hr_point[i];
  	if(((x[j] - LUv[2*j]) < not_zero) || ((LUv[2*j+1] - x[j]) < not_zero)){
  		printf("W[KMATRIX]...\t[KMATRIX_ASL]"
			"Variable \"%d\" (offset %d) has an active bound; sigma = %f.\n", j, i+1, sigma[j]);
			fprintf(somefile, "%d\t%d\t%.g\t%.g\t%.g\t%.g\t%.g\t\t%f\n", j, i+1, z_L[j], z_U[j], LUv[2*j], x[j], LUv[2*j+1], sigma[j]);
  	}
  	positions_rh[j] = i+1;
  	/*printf("j %d, position %d\n", j, positions_rh[j]);*/
  }
  fclose(somefile);


  somefile = fopen("sigma_super_basic.txt", "w");
  for(i=0; i<n_dof; i++){
  	j = hr_point[i];
		fprintf(somefile, "%d\t%d\t%.g\t%.g\t%.g\t%.g\t%.g\t\t%f\n", j, i+1, z_L[j], z_U[j], LUv[2*j], x[j], LUv[2*j+1], sigma[j]); 
  }
  fclose(somefile);

  somefile = fopen("zx.txt", "w");
  for(i=0; i<n_dof; i++){
  	j = hr_point[i];
		fprintf(somefile, "%d\t%d\t%.g\t%.g\n", j, i+1, (LUv[2*j] - x[j]) * z_L[j], (x[j] - LUv[2*j+1]) * z_U[j]); 
  }
  fclose(somefile);

  somefile = fopen("result_red_hess.txt", "w");
  /* fprintf(somefile, "\t%.g", *(x_+ j * K_nrows + hr_point[i])); */
  for(i=0; i<n_dof; i++){
		for(j=0; j<n_dof; j++){
			/*if(j<i){
				fprintf(somefile, "\t%.g", *(x_+ i * K_nrows + hr_point[j]));
			}
			else{
				fprintf(somefile, "\t%.g", *(x_+ j * K_nrows + hr_point[i]));
			}*/
			fprintf(somefile, "\t%.g", (*(x_+ j * K_nrows + hr_point[i]) + *(x_+ i * K_nrows + hr_point[j])) * 0.5 );

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
  if(dot_prod_f != 0){
  	printf("I[KMATRIX]...\t[KMATRIX_ASL]"
		"Dot product preparation phase.\n");
  }
  solve_result_num = 0;
  write_sol(ter_msg, s_star, s_star + n_var, &Oinfo);

  /* evaluate_eigenvalues of the reduced hessian */
  if(eig_rh_eval>0){dsyev_driver(n_dof, x_, K_nrows, hr_point);}
  



  free(positions_rh);
  for(i=0; i<(int)n_r_suff; i++){
  	free(reg_suffix_name[i]);
  }
  free(sigma);
  free(z_L);
	free(z_U);
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
