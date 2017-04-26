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
#include "getstub.h"
#include <assert.h>

#include "sort_by_column.h"
#include "w_append_nz.h"
#include "k_assemble_cc.h"
#include "kmalloc.h"
#include "mc30_driver.h"
//#include "pardiso_driver.h"

#include "get_jac_asl.h"
#include "get_hess_asl.h"

#include "find_inequalities.h"
#include "assemble_corrector_rhs.h"
#include "assemble_rhsds.h"

static int dumm = 1;
static I_Known dumm_kw = {2, &dumm};
static int n_rhs = 0;
static int l_over = 0;
static I_Known l_over_kw = {1, &l_over};

// keywords
static keyword keywds[] = {
  KW("smth", IK_val, &dumm_kw, "Cheers mate the cavalry is here"),
  KW("n_rhs", I_val, &n_rhs, "Number or right hand sides"),	
  KW("no_lambda", IK_val, &l_over_kw, "Override multiplier check(lambda)")
};

static Option_Info Oinfo = 
{"k_matrix", "KKT_matrix_man", "k_matris_opt", keywds, nkeywds};
/*Remember to set reference back*/


int main(int argc, char **argv){
	ASL *asl;
	FILE *f;
	FILE *f_debug;
	/* SufDesc *some_suffix; */
	int i;
	int nnzw; // let's try this
	real *x=NULL, *lambda=NULL;
	char *s=NULL;
	static char sword[] = "rhs_";
	SufDesc *var_f=NULL;
	SufDesc *con_f=NULL;
	SufDesc *rhs_ptr=NULL;
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

	real *b_=NULL, *x_=NULL;
	//real *gf= NULL;
	int *c_flag=NULL;
	char **suf_name=NULL;

	// the pointer to array that will contain
	// the rhs
	real **rhs_baksolve = NULL;

	// The memory allocation
	asl = ASL_alloc(ASL_read_pfgh);


	s = getstops(argv, &Oinfo);

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
//	printf("Number of Right hand sides %d\n", n_rhs);

	// Declare suffix array pointer
	tptr = (SufDecl *)malloc(sizeof(SufDecl)*(n_rhs+2));

	// First two suffixes
	tptr->name = "var_flag";
	tptr->table = 0;
	tptr->kind = ASL_Sufkind_var;
	tptr->nextra = 0;

	(tptr + 1)->name = "con_flag";
	(tptr + 1)->table = 0;
	(tptr + 1)->kind = ASL_Sufkind_con;
	(tptr + 1)->nextra = 0;

	// Suffix for the right hand side creation
	suf_name = (char **)malloc(sizeof(char *)*n_rhs);

	// this bit writes names of the rhs suffixes
	for(i=0; i < n_rhs; i++){
		char numeric_rhs[5];
		numeric_rhs[0] = '\0';
	  sprintf(numeric_rhs, "%d", i);
	  suf_name[i] = (char *)malloc(sizeof(char)*10);
	  *(suf_name[i]) = '\0';
	  strcat(suf_name[i], sword);
	  strcat(suf_name[i], numeric_rhs);
//	  printf("The name of the suffix %s\n", suf_name[i]);
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

	// Declare suffixes
	suf_declare(tptr, (n_rhs + 2));

	// dhis bit setups ASL components e.g. n_var, n_con, etc.
	f = jac0dim(s, 0);

	x = X0 = M1alloc(n_var*sizeof(real));

	// Want initial guess lambda
	want_xpi0 = 2;

	// need to do part of changing sign for y

	pfgh_read(f, ASL_findgroups);
	//printf("%d LBC %f\n",0, Urhsx[2*0]);
	//find_bounds(n_con, LUv);
	
	// change the name of the multipliers to lambda
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
    //return -1;
	}
	if(con_f->u.i == NULL){
		fprintf(stderr, "E[KMATRIX]...\t[KMATRIX_ASL]"
			"u.i empty, no con_flag declared.\n");
	//	return -1;
	}

	//SufDesc *smptr0, *smptr1;

	// smptr0 = suf_get("rhs_0", ASL_Sufkind_con);
	// smptr1 = suf_get("rhs_1", ASL_Sufkind_con);

	for(i=0; i < 1; i++){
	  //printf("The suffix name %s \n", suf_name[i]);
	  rhs_ptr = suf_get(suf_name[i], ASL_Sufkind_con|ASL_Sufkind_real);
	  if((rhs_ptr)->u.r == NULL){
	    fprintf(stderr, "E[KMATRIX]...\t[KMATRIX_ASL]"
	    	"u.r empty, no rhs values declared for rhs_%d.\n", i);
	    return -1;
	  }
        //return -1;
	}
	// }
	for(i=0; i < n_con; i++){
		printf("%d\t%d\n", i, rhs_ptr->u.r[i]);
	}

	c_flag = (int *)malloc(sizeof(int) * n_con);
	
	//constraints flags
	find_ineq_con(n_con, LUrhs, c_flag);

	// Row and colum for the triplet format A matrix
	// size of the number of nonzeroes in the constraint jacobian
	Acol = (fint *)malloc(sizeof(fint)*nzc);
	Arow = (fint *)malloc(sizeof(fint)*nzc);
	Aij  = (real *)malloc(sizeof(real)*nzc);

	fint nerror;
	nerror = 0;
	
	// get_jac_asl function
	// is it better to pass asl as a reference?
	//get_grad_f(asl, x, n_var, n_obj, gf, &nerror);
	get_jac_asl(asl,x,Acol,Arow,Aij,nzc,&nerror);
	get_hess_asl(asl,x,&Wcol,&Wrow,&Wij, n_var, n_con, n_obj, &nnzw, lambda, &nerror);
	//
	printf("I[KMATRIX]...\t[KMATRIX_ASL]"
	"Nonzeroes in the sparse hessian %d\n", nnzw);
	//printf("nonzeroes in the sparse hessian %d\n", nnzw);
	printf("print dummy %d\n", dumm);
	
	if (dumm == 2){
    printf("I[KMATRIX]...\t[KMATRIX_ASL]"
    "Cheers mate, the cavalry is here!");
	}

	printf("size of s %d\n", sizeof(f));
	printf("I[KMATRIX]...\t[KMATRIX_ASL]"
		"Number of Right hand sides %d\n", n_rhs);
	printf("I[KMATRIX]...\t[KMATRIX_ASL]"
		"Number of variables %d\n", n_var);
	printf("I[KMATRIX]...\t[KMATRIX_ASL]"
		"Number of constraints %d\n", n_con);
	
	fint nzA = nzc;
	// Reorder the jacobian
	sortcol(Arow, Acol, Aij, nzA, n_var);
	f_debug = fopen("somefile.txt", "w");
	for(i=0; i < nzA; i++){
    fprintf(f_debug, "%d\t%d\t%.g\n", Arow[i], Acol[i], Aij[i]);
	}
	fclose(f_debug);


	//wnzappnd(Wrow, Wcol, Wij, nnzw, n_var);
	//w_col_sort(Wrow, Wcol, Wij, nnzw, n_var);
	// calculate the space required for the KKT_matrix
	k_space = k_malloc(Wrow, Wcol, Wij, nnzw, n_var, n_con, nzA);

	// temporal space for w_appended
	wa_tsp = k_space - n_con - nzA;
	Krow = (fint *)malloc(sizeof(fint) * k_space);
	Kcol = (fint *)malloc(sizeof(fint) * k_space);
	Kij  = (real *)malloc(sizeof(real) * k_space);
	// wnzappnd(fint *Wrow, fint *Wcol, real *Wij, fint Wnz, fint nvar,
	// fint *Wr_new, fint *Wc_new, real *Wi_new, fint Wnz_new)

	Wr_t = (fint *)malloc(sizeof(fint) * wa_tsp);
	Wc_t = (fint *)malloc(sizeof(fint) * wa_tsp);
	Wi_t = (real *)malloc(sizeof(real) * wa_tsp);

	// append 0's to the main diagonal as necessary
	wnzappnd(Wrow, Wcol, Wij, nnzw, n_var, Wr_t, Wc_t, Wi_t, wa_tsp);
	// fint *Wrow, fint *Wcol, real *Wij, fint Wnz, fint nvar, fint ncon,
	// fint *Arow, fint *Acol, real *Aij, fint Anz
	f_debug = fopen("appended_hessian.txt", "w");
	for(i=0; i<wa_tsp; i++){
    fprintf(f_debug, "\t%ld\t%ld\t%.g\n", Wr_t[i], Wc_t[i], Wi_t[i]);
	}
	fclose(f_debug);


	// the question is if we do not append nz to the W matrix, can k_assemble 
	// still do its work?
	// appearently, yes it will

	// row starts
	K_nrows = n_var + n_con;
	// number of rows of the KKT matrix
	printf("K_nrows %d\n", K_nrows);

	Kr_strt = (fint *)malloc(sizeof(fint) * (K_nrows + 1));
	S_scale = (real *)malloc(sizeof(real) * K_nrows);

	real * Crhs=NULL;
	assemble_corrector_rhs(asl, x, lambda, n_var, n_con, Arow, Acol, Aij, nzA, Crhs, &nerror, c_flag);
	// ASL *asl, real *x, real *lambda,  fint nvar, fint ncon,	fint *Arow, fint *Acol, real *Aij, fint Anz,  real *Crhs, fint *nerror

	k_assemble_cc(Wr_t, Wc_t, Wi_t, wa_tsp, n_var, n_con, 
		Arow, Acol, Aij, nzA, 
		Krow, Kcol, Kij, k_space,
	        Kr_strt);

	mc30driver(K_nrows, k_space, Kij, Krow, Kcol, S_scale);
	//for(i=0; i<K_nrows; i++){
	//        printf("K[%d]=%f\n", i, S_scale[i]);
	//}
	free(S_scale);
	b_ = (real *)malloc(sizeof(real) * K_nrows);
	x_ = (real *)malloc(sizeof(real) * K_nrows);

	for(i=0; i<K_nrows; i++){
	  b_[i] = 1;
	}

	// factorize the matrix
	//	pardiso_driver(Kr_strt, Kcol, Kij, K_nrows, k_space, 1, b_, x_);
	rhs_baksolve = (real **)malloc(sizeof(real *) * n_rhs);
	for(i = 0; i < n_rhs; i++){
		*(rhs_baksolve+i) = (real *)malloc(sizeof(real) * K_nrows);
		//free(*(rhs_baksolve+i));
	}
	assemble_rhsds(n_rhs, 2, rhs_baksolve, n_var, n_con, rhs_ptr);
//int n_rhs, fint rhs_len, 
// real **rhsbksolv, fint nvar, fint ncon, SufDesc *rhs_ptr
	for(i=0; i<n_rhs; i++){
		free(*rhs_baksolve+i);
	}
	free(rhs_baksolve);

	// suf_name = (char **)malloc(sizeof(char *)*n_rhs);

	free(c_flag);
	free(b_);
	free(x_);

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
	//free(rhs_ptr);
	return 0;
}
