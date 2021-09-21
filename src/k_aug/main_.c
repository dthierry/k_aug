///  @file main_.c

/*
 * Created by David Thierry on 4/19/18.
 * @author: David Thierry (dmolinat@andrew.cmu) dav0@lb2016-1
 * Contributors: Flemming Holtorf, Michael Short, Robby Parker, Janghao Park
 * 2/24/2021: @dthierry, I've created some new comments.
*/


#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "getstub.h"
#include "asl.h"


#include "../interfaces/hsl/ma57_driver.h"
#include "get_jac_asl_aug.h"
#include "get_hess_asl_aug.h"
#include "get_grad_f.h"
#include "find_inequalities.h"
#include "assemble_rhs_rh.h"
#include "assemble_rhs_dcdp.h"
#include "suffix_decl_hand.h"
#include "csr_driver.h"
#include "sigma_compute.h"
#include "mu_adjust_primal.h"
#include "../matrix/dsyev_driver.h"
#include "../matrix/dpotri_driver.h"


#include "k_aug_data.h"
#include "../../config_kaug.h"
#ifdef USE_MC30
#include "../HSL/mc30_driver.h"
#else
#include "../interfaces/hsl/mc19_driver.h"
#endif

#define NUM_REG_SUF 8


static double not_zero = 1.84e-04;  ///< Parameter used to print warning about active bounds (initial value).

static int n_rhs = 0;   ///< We used to have to manually set the number of rhs, not anymore. This is useless now.
static real log10mu = -11.0;

/* Available keyword options
 * Option Description */
static char sc_comp_inv[] = {"compute_inv"}; /* The actual Reduced Hessian */
static char sc_comp_inv_verb[] = {"Compute the inv(inv(red_hess)) i.e. the Red. Hess"};

static char sc_dbg_kkt[] = {"deb_kkt"}; /* Debugging */
static char sc_dbg_kkt_verb[] = {"Override dof checking, for debugging purposes"};

static char sc_dot_pr_f[] = {"dot_prod"}; /* No use now. */
static char sc_dsdpmode[] = {"dsdp_mode"};  /* Like in sipopt. */
static char sc_dsdp_mode_nverb[] = {"Compute the dsdp for constraints of kind C(x) - P = 0 (linear P)"};

static char name1[] = {"smth"};  /* Nothing here */

static char sc_e_eval[] = {"eig_rh"};
static char sc_e_eval_verb[] = {"Compute the eigenvalues of the inv of red_hessian"};

static char sc_n_rhsopt_[] = {"n_rhs"};  /* No use now. */
static char sc_no_barrieropt_[] = {"no_barrier"};  /* Skips adding sigma_i in the Hessian block of the KKT */
static char sc_no_lambdaopt_[] = {"no_lambda"};
static char sc_no_scaleopt_[]  = {"no_scale"};  /* Skips the scaling of the KKT */
static char sc_not_zero[] = {"not_zero"};  /* Use this to change some warnings */

static char sc_print_kkt[] = {"print_kkt"};  /* Useful if you want the kkt matrix and so forth. */
static char sc_print_kkt_noverb[] = {"Prints the contents of the kkt matrix in a file."};

static char sc_square_override[] = {"square_problem"};  /* Do not compute sigmas */
static char sc_square_override_noverb[] = {"If nvar = neq, assume there are no bounds and sigma = 0"};

static char sc_no_inertia[] = {"no_inertia"};
static char sc_no_inertia_verb[] = {"Skip inertia correction"};

static char sc_target_log10mu[] = {"target_log10mu"};
static char sc_target_log10mu_verb[] = {"Target value for log10mu"};

/*
 * Initialization of some numeric types of optional values.
 */

static int dumm = 1;
static I_Known dumm_kw = {2, &dumm};

static int l_over = 0;
static I_Known l_over_kw = {1, &l_over};

static int deb_kkt = 0;
static I_Known deb_kkt_kw = {1, &deb_kkt};

static int dot_prod_f = 0;
static I_Known dot_p_kw = {1, &dot_prod_f};

static int eig_rh_eval = 0;
static I_Known e_eval_kw = {1, &eig_rh_eval};

static int comp_inv_ = 0;
static I_Known comp_inv_kw = {1, &comp_inv_};

static int no_barrier = 1;
static I_Known nbarrier_kw = {0, &no_barrier};

static int no_scale = 1;
static I_Known nscale_kw = {0, &no_scale};

static int dsdp_mode = 0;
static I_Known dsdp_mode_kw = {1, &dsdp_mode};

static int no_inertia = 0;
static I_Known no_inertia_kw = {1, &no_inertia};

static int print_kkt = 0;
static I_Known print_kkt_kw = {1, &print_kkt};

static int square_override = 0;
static I_Known square_override_kw = {1, &square_override};

/*static char dof_v[] = {"dof_v"};*/

/*
 * keywords, they must be in alphabetical order!
 * */
static keyword keywds[] = {
        KW(sc_comp_inv, IK_val, &comp_inv_kw, sc_comp_inv_verb),
        KW(sc_dbg_kkt, IK_val, &deb_kkt_kw, sc_dbg_kkt_verb),
        KW(sc_dot_pr_f, IK_val, &dot_p_kw, sc_dot_pr_f),
        KW(sc_dsdpmode, IK_val, &dsdp_mode_kw, sc_dsdp_mode_nverb),
        KW(sc_e_eval, IK_val, &e_eval_kw, sc_e_eval_verb),
        KW(sc_n_rhsopt_ , I_val, &n_rhs, sc_n_rhsopt_),
        KW(sc_no_barrieropt_ , IK_val, &nbarrier_kw, sc_no_barrieropt_),
        KW(sc_no_inertia, IK_val, &no_inertia_kw, sc_no_inertia_verb),
        KW(sc_no_lambdaopt_ , IK_val, &l_over_kw, sc_no_lambdaopt_),
        KW(sc_no_scaleopt_ , IK_val, &nscale_kw, sc_no_scaleopt_),
        KW(sc_not_zero , D_val, &not_zero, sc_not_zero),
        KW(sc_print_kkt, IK_val, &print_kkt_kw, sc_print_kkt_noverb),
        KW(name1 , IK_val, &dumm_kw, name1),  /* This does not do anything */
        KW(sc_square_override, IK_val, &square_override_kw, sc_square_override_noverb),
        KW(sc_target_log10mu , D_val, &log10mu, sc_target_log10mu_verb),
};


static char banner[] = {"[K_AUG] 0.1.0, Part of the IDAES PSE framework\n"
                        "Please visit https://idaes.org/\n"};
static char sc_k_[] = {"k_aug"};
static char sc_k_o_[] = {"k_aug_options"};

static Option_Info Oinfo;



/**
 * Main function, Reads nl file, allocates data structures, calls assembling funcs.
 * @param argc Typically the .nl file followed by the options.
 * @param argv
 * @return
 */
int main(int argc, char **argv){
    ASL *asl;
    FILE *f;

    /*
     * Sometimes we use ASL's types, sometimes native types.
     * */
    int i, j, k;
    int n_dof=0;
    int n_vx=0; /* variable of interest for dxdp*/
    fint nnzw; /* let's try this */
    real *x=NULL, *lambda=NULL;
    char *s=NULL;
    /*
     * ASL's SufDesc *some_suffix;
     */
    SufDesc *var_f=NULL;
    SufDesc *suf_zL = NULL;
    SufDesc	*suf_zU = NULL;
    SufDesc *f_timestamp = NULL;
    SufDesc *dcdp = NULL;
    SufDesc *var_order_suf = NULL;
    SufDesc *con_order_suf = NULL;
    SufDesc **rhs_ptr=NULL;
    SufDecl *suf_ptr=NULL;

    /*
     * These are all ASL's types. Hopefully ASL is always configured correctly.
     */
    real *z_L=NULL, *z_U=NULL, *sigma=NULL;
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
    real *gf = NULL;
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

    char ter_msg[] = {"I[K_AUG]...[K_AUG_ASL]"
                      "All done."};

    unsigned n_r_suff = NUM_REG_SUF;

    /*
     * Suffixes names. Add new suffixes names here.
     * I decided to just use a sequence of _sfx_X for no good reason.
     */
    char _sfx_1[] = {"dof_v"};
    char _sfx_2[] = {"rh_name"};
    char _sfx_3[] = {"ipopt_zL_in"};
    char _sfx_4[] = {"ipopt_zU_in"};
    char _sfx_5[] = {"f_timestamp"};
    char _sfx_6[] = {"dcdp"};
    char _sfx_7[] = {"var_order"};
    char _sfx_8[] = {"con_order"};

    int _is_not_irh_pdf=1;
    clock_t start_c, ev_as_kkt_c, fact_kkt_c, total_c;
    double ev_as_time, fact_time, overall_time;
    time_t timestamp;
    char _chr_timest[15] = ""; /* The timestamp */
    char _file_name_[30] = ""; /* */
    char WantModifiedJac = 0;
    int *vModJac = NULL, *cModJac = NULL;
    double logmu0;
    struct stat st;

    /*
     * Inertia data-structures
     */
    inertia_params inrt_parms;
    inertia_options inrt_opts;
    inertia_perts inrt_pert;
    linsol_opts ls_opts;
    /*
     * Instability made me do this.
     */
    putenv("OMP_NUM_THREADS=1");
    /*
     * Inertia parameters, Inspired by ipopt's parameters.
     */
    inrt_parms.dmin = 10e-20;
    inrt_parms.dmax = 10e+40;
    inrt_parms.d_w0 = 1e-08;
    inrt_parms.kbp = 100;
    inrt_parms.kp = 8;
    inrt_parms.km = 1.0/3.0;

    inrt_parms.dcb = 1e-08;
    inrt_parms.kc = 0.25;


    timestamp = time(NULL);
    start_c = clock();

    inrt_opts.always_perturb_jacobian = 0;
    inrt_opts.no_inertia = no_inertia;

    inrt_pert.d_c = 0.0;
    inrt_pert.d_w = 0.0;
    inrt_pert.d_c_last = 0.0;
    inrt_pert.d_w_last = 0.0;
    inrt_pert.jacobian_perturbed = 0;

    ls_opts.want_accurate = 1;
    ls_opts.max_inertia_steps = 20;
    ls_opts.max_refinement_steps = 5;
    ls_opts.residual_ratio_max = 1e-10;
    ls_opts.max_memory_al_steps = 10;
    ls_opts.pivot_tol0 = 1e-06;
    ls_opts.pivtol_max = 0.5;

    /*
     * ASL stuff.
     */
    Oinfo.sname = sc_k_;
    Oinfo.bsname = banner;
    Oinfo.opname = sc_k_o_;
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

    /*
     * In case an user wants to get the debug files
     */
#ifndef PRINT_VERBOSE
    if (stat("./kaug_debug", &st) == -1){
        printf("I[K_AUG]...\t[K_AUG_ASL]Creating kaug_debug folder\n");
#ifndef USE_MINGW64
        mkdir("./kaug_debug", 0700);
#else
        mkdir("./kaug_debug");
        /* David: For some reason MINGW uses this signature */
#endif
    }
#endif
    /* The memory allocation for asl data structure */
    asl = ASL_alloc(ASL_read_pfgh);
    /*
     * THE Read.
     */
    s = getstops(argv, &Oinfo);


    if (!s) {
        printf("W[K_AUG]...\t[K_AUG_ASL]"
               "No input\n");
        return 1;
    }
    else {
        printf("I[K_AUG]...\t[K_AUG_ASL]"
               "File read successfully.\n");
    }

    if (n_rhs == 0){
        fprintf(stderr, "W[K_AUG]...\t[K_AUG_ASL]"
                        "No n_rhs declared.\n");
    }

    if (no_inertia){
        fprintf(stderr, "W[K_AUG]...\t[K_AUG_ASL]"
                        "Inertia correction skip.\n");
    }

    fprintf(stderr, "W[K_AUG]...\t[K_AUG_ASL]"
                    "Target log10mu:= %.g.\n", log10mu);

    if (l_over){
        fprintf(stderr, "W[K_AUG]...\t[K_AUG_ASL]"
                        "Multiplier check override.\n");
    }
    if (dot_prod_f){
        fprintf(stderr, "W[K_AUG]...\t[K_AUG_ASL]"
                        "Dot product preparation.\n");
    }


    /*
     * Allocate suffix names (regular)
     */
    reg_suffix_name =   (char **)malloc(sizeof(char *) * n_r_suff);
    for(i=0; i < (int)n_r_suff; i++){
        reg_suffix_name[i] = (char *)malloc(sizeof(char) * 16 );
        reg_suffix_name[i][0] = '\0';
    }

    strcat(reg_suffix_name[0], _sfx_1);
    strcat(reg_suffix_name[1], _sfx_2);
    strcat(reg_suffix_name[2], _sfx_3);
    strcat(reg_suffix_name[3], _sfx_4);
    strcat(reg_suffix_name[4], _sfx_5);
    strcat(reg_suffix_name[5], _sfx_6);
    strcat(reg_suffix_name[6], _sfx_7);
    strcat(reg_suffix_name[7], _sfx_8);

    if(n_rhs > 0){
        suf_ptr = (SufDecl *)malloc(sizeof(SufDecl)*(n_rhs + n_r_suff));
        rhs_name = (char **)malloc(sizeof(char *)*n_rhs);
        for(i=0; i<n_rhs; i++){
            rhs_name[i] = (char *)malloc(sizeof(char) * 32); /* 32 bit long digit;
			 */
        }
        /*
         * Declaration of suffixes
         */
        suffix_decl_hand(suf_ptr, reg_suffix_name, rhs_name, n_r_suff, n_rhs);
    }
    else{
        suf_ptr = (SufDecl *)malloc(sizeof(SufDecl) * n_r_suff);
        suffix_decl_hand(suf_ptr, reg_suffix_name, rhs_name, n_r_suff, n_rhs);
    }


    /* Declare suffixes */
    if(n_rhs > 0){suf_declare(suf_ptr, (n_rhs + n_r_suff));}
    else{suf_declare(suf_ptr, n_r_suff);}


    /*
     * This bit setups ASL components e.g. n_var, n_con, etc.
     */
    f = jac0dim(s, (fint)strlen(s));


    printf("I[K_AUG]...\t[K_AUG_ASL]"
           "Number of Right hand sides: %d\n", n_rhs);
    printf("I[K_AUG]...\t[K_AUG_ASL]"
           "Number of variables       : %d\n", n_var);
    printf("I[K_AUG]...\t[K_AUG_ASL]"
           "Number of constraints     : %d\n", n_con);
    printf("I[K_AUG]...\t[K_AUG_ASL]"
           "Number of valid n_dof     : %d\n", n_var - n_con );


    if ((n_var - n_con) < 0){
        if(deb_kkt>0){
            fprintf(stderr, "W[K_AUG]...\t[K_AUG_ASL]"
                            "KKT check!\n");
        }
        else{
            printf("E[K_AUG]...\t[K_AUG_ASL]"
                   "nvar < ncon. This problem is not valid.\n");
            exit(-1);
        }
    }
    else if (n_var == n_con){
        fprintf(stderr, "W[K_AUG]...\t[K_AUG_ASL]"
                        "This problem has no degrees of freedom\n"
                        "Pass the option square_override for the desired behaviour\n");
        if(deb_kkt>0){
            fprintf(stderr, "W[K_AUG]...\t[K_AUG_ASL]"
                            "KKT check!\n");
        }
    }
    /*
     * Primal-Dual Values
     */
    x 		 = X0  = M1alloc(n_var*sizeof(real));
    lambda = pi0 = M1alloc(n_con*sizeof(real));
    obj_no = 0;


    pfgh_read(f, ASL_findgroups);

    /*
     * In case we don't declare the multipliers.
     * Would this create a memory leak? This is pretty old section.
     */
    if(lambda==NULL && l_over == 0){
        printf("E[K_AUG]...\t[K_AUG_ASL]"
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

    /*
     * Get the values of suffixes.
     */
    var_f = suf_get(reg_suffix_name[0], ASL_Sufkind_var);

    suf_zL = suf_get(reg_suffix_name[2], ASL_Sufkind_var| ASL_Sufkind_real);
    suf_zU = suf_get(reg_suffix_name[3], ASL_Sufkind_var| ASL_Sufkind_real);

    f_timestamp = suf_get(reg_suffix_name[4], ASL_Sufkind_prob);
    dcdp = suf_get(reg_suffix_name[5], ASL_Sufkind_con);

    var_order_suf = suf_get(reg_suffix_name[6], ASL_Sufkind_var);
    con_order_suf = suf_get(reg_suffix_name[7], ASL_Sufkind_con);

    z_L = (real *)malloc(sizeof(real) * n_var);
    z_U = (real *)malloc(sizeof(real) * n_var);


    memset(z_L, 0, sizeof(real) * n_var);
    memset(z_U, 0, sizeof(real) * n_var);


    if(!(suf_zL->u.r)){
        fprintf(stderr, "W[K_AUG_ASL]...\t[K_AUG_ASL]"
                        "No ipopt_zL_in suffix declared, setting zL = 0.\n");
    }
    else{
        for(i=0; i< n_var; i++){
            z_L[i] = suf_zL->u.r[i];
        }
    }
    if(!(suf_zU->u.r)){
        fprintf(stderr, "W[K_AUG_ASL]...\t[K_AUG_ASL]"
                        "No ipopt_zU_in suffix declared, setting zU = 0.\n");
    }
    else{
        for(i=0; i< n_var; i++){
            z_U[i] = suf_zU->u.r[i];
        }
    }


    if (dsdp_mode > 0) {
        strcat(_file_name_, "dsdp_in_");
    }
    else{
        strcat(_file_name_, "dot_in_");
    }
    if(!(f_timestamp->u.i)){
        fprintf(stderr, "W[K_AUG]...\t[K_AUG_ASL]"
                        "No f_timestamp suffix declared, Fallback to default writing mode.\n");
    }
    else{
        printf("I[K_AUG]...\t[K_AUG_ASL]"
               "Timestamp suffix = %d.\n\n", *(f_timestamp->u.i));
        sprintf(_chr_timest, "%d", *(f_timestamp->u.i));
        /*fprintf(stderr, "This goes here %s\n", _chr_timest);*/
    }
    strcat(_file_name_, _chr_timest);
    strcat(_file_name_, ".in");
    fprintf(stderr, "I[K_AUG]...\t[K_AUG_ASL] Filename for dot_sens %s\n", _file_name_);

    /*
     * This directives print some debugging files.
     */

#ifndef PRINT_VERBOSE
    somefile = fopen("./kaug_debug/primal0.txt", "w");
    for(i=0; i< n_var; i++){
        fprintf(somefile, "%.g\n", x[i]);
    }
    fclose(somefile);
#endif
    /*
     * Use the mu value to correct the slack/multipliers.
     * This should be connected somehow to Ipopt, (which I have tested through suffixes).
     * We need to incorporate those changes to the ipopt/asl interface.
     */
    logmu0 = log10mu;
    if (!square_override){
        mu_adjust_x(n_var, x, LUv, z_L, z_U, log10mu, &logmu0);
    }

#ifndef PRINT_VERBOSE
    somefile = fopen("./kaug_debug/primal1.txt", "w");
    for(i=0; i< n_var; i++){fprintf(somefile, "%.g\n", x[i]);}
    fclose(somefile);
#endif
    /*
     * sigma_i = z_i/x_i
     */
    sigma = (real *)malloc(sizeof(real) * n_var);
    memset(sigma, 0, sizeof(real) * n_var);
    if (print_kkt){
        printf("yes I hear you %d\n", print_kkt);
    }
    /* Check if we do red_hess, deb_kkt or dsdp*/
    if (deb_kkt > 0) {
        fprintf(stderr, "W[K_AUG]...\t[K_AUG_ASL]"
                        "KKT check!\n");
    } else if (dsdp_mode > 0) {
        fprintf(stderr, "W[K_AUG]...\t[K_AUG_ASL]"
                        "dsdp for linear C(x) + I*p = 0 override.\n");
        if (dcdp->u.r == NULL && dcdp->u.i == NULL) {
            fprintf(stderr, "E[K_AUG]...\t[K_AUG_ASL]"
                            "suffix empty no dcdp declared!\n");
            exit(-1);
        }
        somefile = fopen("conorder.txt", "w");
        for (i = 0; i < n_con; i++) { fprintf(somefile, "%d\n", dcdp->u.i[i]); }
        fclose(somefile);
    } else if (print_kkt) {
        printf("W[k_aug]...\t[k_aug_asl]"
               "printing kkt matrix.\n");
    } else if (var_f->u.r == NULL && var_f->u.i == NULL) {
        printf("W[k_aug]...\t[k_aug_asl]"
               "DEFAULT mode!.\n");
        fprintf(stderr, "E[K_AUG]...\t[K_AUG_ASL]"
                        "suffix empty no n_dof declared!\n");
        exit(-1);
    }

    /*
     * Actually compute sigma
     */
    if(!square_override) {
        compute_sigma(asl, n_var, x, z_L, z_U, sigma, logmu0);
    }


    /* todo: This will leak memory. */
    if(n_rhs > 1){
        rhs_ptr = (SufDesc **)malloc(sizeof(SufDesc *) * n_rhs);
        for(i=0; i < n_rhs; i++){
            *(rhs_ptr + i)= suf_get(rhs_name[i], ASL_Sufkind_con|ASL_Sufkind_real);
            if((*(rhs_ptr + i))->u.r == NULL){
                fprintf(stderr, "E[K_AUG]...\t[K_AUG_ASL]"
                                "No rhs values declared for rhs_%d.\n", i);
                exit(-1);
            }
        }
    }

    c_flag = (int *)malloc(sizeof(int) * n_con); /* Flags for ineq or equalities*/

    /*
     * Constraints flags
     */
    find_ineq_con(n_con, LUrhs, c_flag); /* Find the inequality constraints */

    /* Row and column for the triplet format A matrix */
    /* size of the number of nonzeroes in the constraint jacobian */
    Acol = (fint *)malloc(sizeof(fint)*nzc);
    Arow = (fint *)malloc(sizeof(fint)*nzc);
    Aij  = (real *)malloc(sizeof(real)*nzc);


    nerror = 0;
    printf("I[K_AUG]...\t[K_AUG_ASL]"
           "Nonzeroes in the sparse Jacobian %d\n", nzc);


    /*
     * Ask for the J and H values
     */
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
        printf("I[K_AUG]...\t[K_AUG_ASL]"
               "Barrier term added.\n");
    }
    /*
     * Sometimes we wanted to do wacky stuff with the matrices.
     * E.g., get a permuted Jacobian or even partition it.
     * This was where we were doing this.
     * Not relevant for most users anyway.
     */
    if(deb_kkt > 0){
        if(var_order_suf->u.i != NULL && con_order_suf->u.i != NULL){
            WantModifiedJac = 1;
            vModJac = var_order_suf->u.i;
            cModJac = con_order_suf->u.i;
            printf("I[K_AUG]...\t[K_AUG_ASL]"
                   "var_order & con_order suffixes found.\n");
            somefile = fopen("orders_v_.txt", "w");
            for(i=0; i<n_var; i++){fprintf(somefile, "%d\n",*(vModJac + i));}
            fclose(somefile);
            somefile = fopen("orders_c_.txt", "w");
            for(i=0; i<n_con; i++){fprintf(somefile, "%d\n",*(cModJac + i));}
            fclose(somefile);
        }
        else if(var_order_suf->u.i != NULL){
            WantModifiedJac = 2;
            vModJac = var_order_suf->u.i;
            printf("I[K_AUG]...\t[K_AUG_ASL]"
                   "var_order suffix found.\n");
            somefile = fopen("orders_v_.txt", "w");
            for(i=0; i<n_var; i++){fprintf(somefile, "%d\n",*(vModJac + i));}
            fclose(somefile);
        }
        else if(con_order_suf->u.i != NULL){
            WantModifiedJac = 3;
            cModJac = con_order_suf->u.i;
            printf("I[K_AUG]...\t[K_AUG_ASL]"
                   "con_order suffix found.\n");
            somefile = fopen("orders_c_.txt", "w");
            for(i=0; i<n_con; i++){fprintf(somefile, "%d\n",*(cModJac + i));}
            fclose(somefile);
        }

        if(WantModifiedJac > 0){
            printf("I[K_AUG]...\t[K_AUG_ASL]"
                   "Write modified Jacobian file.\n");
            printf("I[K_AUG]...\t[K_AUG_ASL]"
                   "Please verify your suffixes with the orders_x.txt files.\n");
            somefile = fopen("testfile.dat", "w");
            for(i=0; i<nzc; i++){fprintf(somefile, "%d\t%d\t%.g\n", Acol[i], Arow[i], Aij[i]);}

            fclose(somefile);
            somefile = fopen("modified_jacobian.dat", "w");
            /* Jacobian starts at 1 [MATLAB] */
            /* The downside is that if you make a mistake on the suffix values, you won't know if the
               modified jacobian is correct. */
            if(WantModifiedJac == 1){
                for(i=0; i<nzc; i++){
                    fprintf(somefile, "%d\t%d\t%.g\n", *(cModJac+Acol[i]-1), *(vModJac+Arow[i]-1), Aij[i]);
                }
            }
            else if(WantModifiedJac == 2){
                for(i=0; i<nzc; i++){
                    fprintf(somefile, "%d\t%d\t%.g\n", Acol[i]+1, *(vModJac+Arow[i]-1), Aij[i]);
                }
            }
            else if(WantModifiedJac == 3){
                for(i=0; i<nzc; i++){
                    fprintf(somefile, "%d\t%d\t%.g\n", *(cModJac+Acol[i]-1), Arow[i]+1, Aij[i]);
                }
            }
            fclose(somefile);
        }

        solve_result_num = 0;
        write_sol(ter_msg, x, lambda, &Oinfo);
        ASL_free(&asl);
        free(c_flag);
        free(z_L);
        free(z_U);
        free(sigma);
        free(Acol);
        free(Arow);
        free(Aij);
        free(Wcol);
        free(Wrow);
        free(Wij);
        free(nz_row_a);
        free(nz_row_w);
        free(md_off_w);
        for(i=0; i<(int)n_r_suff; i++){
            free(reg_suffix_name[i]);
        }
        free(reg_suffix_name);
        free(suf_ptr);
        /* Execution should end here, no memory leaks */
        return 0;}
    nzA = nzc;
    nzW = nnzw + miss_nz_w;
    /* Recomputed number of nz in the Hessian-of-Lagrangian */

    /*
     * Generate the data structure for K
     */
    csr_driver(n_var, n_con, nzW, nzA, nz_row_w, nz_row_a,
               Wrow, Wcol, Wij, Arow, Acol, Aij,
               &Krow, &Kcol, &Kij, &Kr_strt);


    K_nrows = n_var + n_con; /* Number of rows of the KKT matrix (no ineq) */
    nzK = nzA + nzW + n_con; /* NZ in the kkt matrix (for pardiso, yes)*/
    assert(Krow != NULL);
    ev_as_kkt_c = clock();

    /*
     * I got asked to do this. Print the matrices to files.
     */
    if(print_kkt){
        if (stat("./GJH", &st) == -1){
#ifndef USE_MINGW64
            mkdir("./GJH", 0700); /* Some versions of MINGW have a different signature :S */
#else
            mkdir("./GJH");
#endif
            printf("I[K_AUG]...\t[K_AUG_ASL]"
                   "New folder created ./GJH !...\n");
        }
        gf = (real *)calloc(sizeof(real), n_var);
        get_grad_f(asl, x, n_var, gf, &nerror);

        somefile = fopen("./GJH/gradient_f_print.txt", "w");
        for(k=0; k<n_var; k++){fprintf(somefile, "%.g\n", gf[k]);}
        fclose(somefile);

        somefile = fopen("./GJH/A_print.txt", "w");
        for(k=0; k<nzA; k++){fprintf(somefile, "%d\t%d\t%.g\n", Arow[k], Acol[k], Aij[k]);}
        fclose(somefile);

        somefile = fopen("./GJH/H_print.txt", "w");
        for(k=0; k<nzW; k++){fprintf(somefile, "%d\t%d\t%.g\n", Wrow[k], Wcol[k], Wij[k]);}
        fclose(somefile);

        somefile = fopen("./GJH/kkt_print.txt", "w");
        for(k=0; k<nzK; k++){fprintf(somefile, "%d\t%d\t%.g\n", Krow[k], Kcol[k], Kij[k]);}
        fclose(somefile);
        solve_result_num = 0;
        write_sol(ter_msg, x, lambda, &Oinfo);

        ASL_free(&asl);
        free(gf);
        free(c_flag);
        free(z_L);
        free(z_U);
        free(sigma);
        free(Acol);
        free(Arow);
        free(Aij);
        free(Wcol);
        free(Wrow);
        free(Wij);
        free(nz_row_a);
        free(nz_row_w);
        free(md_off_w);
        free(Krow);
        free(Kcol);
        free(Kij);
        free(Kr_strt);

        for(i=0; i<(int)n_r_suff; i++){
            free(reg_suffix_name[i]);
        }
        free(reg_suffix_name);
        free(suf_ptr);
        /*
         * No memory leaks.
         */
        return 0;
    }


    S_scale = (real *)calloc(sizeof(real), K_nrows);

    /*
     * MC30 should not even be used anymore.
     */
#ifdef USE_MC30
    printf("I[K_AUG]...\t[K_AUG_ASL]"
           "MC30 scaling...\n");
    mc30driver(K_nrows, nzK, Kij, Krow, Kcol, S_scale);
#else
    printf("I[K_AUG]...\t[K_AUG_ASL]"
           "MC19 scaling...\n");
    mc19driver(K_nrows, nzK, Kij, Krow, Kcol, S_scale);
#endif


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

#ifndef PRINT_VERBOSE
    somefile = fopen("./kaug_debug/primal_dual.txt", "w");
    for(i=0; i<K_nrows; i++){
        fprintf(somefile, "\t%f\n", s_star[i]);
    }
    fclose(somefile);
#endif

    /*
     * This is where we assemble the rhs. Note that There are no derivatives in the RHS only 1's and 0's
     */
    /* problem: all stuff associated with n_dof, solution: use it again for dsdp*/
    if (dsdp_mode > 0) {
        assemble_rhs_dcdp(&rhs_baksolve, n_var, n_con, &n_dof, &n_vx, dcdp, &hr_point, var_order_suf);
    }
    else{
        assemble_rhs_rh(&rhs_baksolve, n_var, n_con, &n_dof, var_f, &hr_point);
    }

    x_           = (real *)calloc(K_nrows * (n_dof), sizeof(real));
    positions_rh = (int *)malloc(n_var * sizeof(int));


    /*
     * Scale matrix & rhs
     */
    if(no_scale > 0){
//#ifdef USE_MC30
        for(i=0; i< nzK; i++){
            Kij[i] = Kij[i] * exp(S_scale[Kcol[i]-1]) * exp(S_scale[Krow[i]-1]);
            /* Kij[i] = Kij[i] * exp((S_scale[Krow[i]-1] + C_scale[Kcol[i]-1])/2.); */
        }
        for(i=0; i< n_dof; i++){
            for(j=0; j < K_nrows; j++){
                *(rhs_baksolve + i*K_nrows + j) =
                        *(rhs_baksolve + i*K_nrows + j) * exp(S_scale[j]);
            }
        }
#ifndef PRINT_VERBOSE
        f = fopen("./kaug_debug/scale_facts.txt", "w");
        somefile = fopen("./kaug_debug/rhs_sens_scaled", "w");
        for(j=0; j < K_nrows; j++){
            for(i=0; i < n_dof; i++){
                fprintf(f, "%.g\n", S_scale[j]);
                fprintf(somefile, "\t%f\t", *(rhs_baksolve + i*K_nrows + j) );
            }
            fprintf(somefile, "\n");
        }
        fclose(somefile);
        fclose(f);
#endif
    }
    else{
        fprintf(stderr, "W[K_AUG]...\t[K_AUG_ASL]"
                        "The scaling has been skipped. \n");
    }

    /* Factorize the matrix
     * I want to dynamically load different solvers, there should be a generic solver interface.
     */
    /* mumps_driver(Kr_strt, Krow, Kcol, Kij, K_nrows, n_dof, rhs_baksolve, x_, n_var, n_con, no_inertia, nzK,
                 &inrt_pert, inrt_parms, &inrt_opts, logmu0, ls_opts); */
    ma57_driver(Kr_strt, Krow, Kcol, Kij, K_nrows, n_dof, rhs_baksolve, x_, n_var, n_con, no_inertia, nzK,
                &inrt_pert, inrt_parms, &inrt_opts, logmu0, ls_opts);

    /*x_ = rhs_baksolve;*/
    printf("I[K_AUG]...\t[K_AUG_ASL]"
           "Linear solver done. \n");
    fact_kkt_c = clock();

    /*
     * Only for debugging.
     */
#ifndef PRINT_VERBOSE
    somefile = fopen("./kaug_debug/result_lin_sol.txt", "w");
    for(i=0; i<K_nrows; i++){
        fprintf(somefile, "\t%d", i);
        for(j=0; j<n_dof; j++){
            fprintf(somefile, "\t%f", *(x_+ j * K_nrows + i));
        }
        fprintf(somefile, "\n");
    }
    fclose(somefile);
#endif
#ifndef PRINT_VERBOSE
    somefile = fopen("./kaug_debug/result_unscaled.txt", "w");
#endif
    if(no_scale > 0){
//#ifdef USE_MC30
        for(i=0; i<K_nrows; i++){
            for(j=0; j<n_dof; j++){
                *(x_+ j * K_nrows + i) = *(x_+ j * K_nrows + i) * exp(S_scale[i]);
#ifndef PRINT_VERBOSE
                fprintf(somefile, "\t%f", *(x_+ j * K_nrows + i));
#endif

            }
#ifndef PRINT_VERBOSE
            fprintf(somefile, "\n");
#endif
        }
//#endif
    }
#ifndef PRINT_VERBOSE
    fclose(somefile);
#endif

    if(var_order_suf->u.i){
        printf("I[K_AUG]...\t[K_AUG_ASL]"
               "var_order suffix detected. Make sure var_order[i] > 0.\n");
#ifndef PRINT_VERBOSE
        somefile = fopen("./kaug_debug/varorder.txt", "w");
        for(i=0; i<n_var; i++){fprintf(somefile, "%d\n", *(var_order_suf->u.i + i));}
        fclose(somefile);
#endif
        somefile = fopen("./dxdp_.dat", "w");

        for(i=0; i<n_vx; i++){
            j = hr_point[i]; /* The row */
            /*fprintf(somefile, "%d", j);*/
            for(k=0; k<n_dof; k++){fprintf(somefile, "\t%.g", *(x_+ k * K_nrows + j));}
            fprintf(somefile, "\n");
        }

        fclose(somefile);
    }
    /* This should write the dsdp file as well. */
    somefile = fopen(_file_name_, "w"); /* For dot_driver */
    for(i=0; i<n_dof; i++){
        for(j=0; j<K_nrows; j++){
            fprintf(somefile, "\t%.g\n", *(x_+ i * K_nrows + j));
        }
    }
    fclose(somefile);
    /* Skip this if dsdp */
    if (dsdp_mode > 0) {}
    else{
        memset(positions_rh, 0, sizeof(int)*n_var);
        somefile = fopen("sigma_warnings.txt", "w");
        for(i=0; i<n_dof; i++){
            j = hr_point[i];
            if(((x[j] - LUv[2*j]) < not_zero) || ((LUv[2*j+1] - x[j]) < not_zero)){
                fprintf(stderr, "W[K_AUG]...\t[K_AUG_ASL]"
                                "Variable \"%d\" (offset %d) has an active bound; sigma = %f.\n",
                        j, i+1, sigma[j]);
                fprintf(somefile, "%d\t%d\t%.g\t%.g\t%.g\t%.g\t%.g\t\t%f\n",
                        j, i+1, z_L[j], z_U[j], LUv[2*j], x[j], LUv[2*j+1], sigma[j]);
            }
            positions_rh[j] = i;
            /*printf("j %d, position %d\n", j, positions_rh[j]);*/
        }
        fclose(somefile);
#ifndef PRINT_VERBOSE
        somefile = fopen("./kaug_debug/sigma_super_basic.txt", "w");
        for(i=0; i<n_dof; i++){
            j = hr_point[i];
            fprintf(somefile, "%d\t%d\t%.g\t%.g\t%.g\t%.g\t%.g\t\t%f\n",
                    j, i+1, z_L[j], z_U[j], LUv[2*j], x[j], LUv[2*j+1], sigma[j]);
        }
        fclose(somefile);

        somefile = fopen("./kaug_debug/zx.txt", "w");
        for(i=0; i<n_dof; i++){
            j = hr_point[i];
            fprintf(somefile, "%d\t%d\t%.g\t%.g\n",
                    j, i+1, (LUv[2*j] - x[j]) * z_L[j], (x[j] - LUv[2*j+1]) * z_U[j]);
        }
        fclose(somefile);
#endif
        /*
         * Get some file with results
         */
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
                fprintf(somefile, "\t%.g",
                        (*(x_+ j * K_nrows + hr_point[i]) + *(x_+ i * K_nrows + hr_point[j])) * 0.5 );

            }
            fprintf(somefile, "\n");
        }
        fclose(somefile);
    }



#ifndef PRINT_VERBOSE
    somefile = fopen("./kaug_debug/result_primal_dual.txt", "w");
    for(i=0; i<K_nrows; i++){
        fprintf(somefile, "\t%f\n", s_star[i]);
    }
    fclose(somefile);
#endif

    /*
     * Write the permutation of the RedHess
     */
    suf_iput(reg_suffix_name[1], ASL_Sufkind_var, positions_rh);

    if(dot_prod_f != 0){
        printf("I[K_AUG]...\t[K_AUG_ASL]"
               "Dot product preparation phase.\n");
    }
    solve_result_num = 0;

    /*
     * Write back
     */
    write_sol(ter_msg, s_star, s_star + n_var, &Oinfo);


    /* evaluate_eigenvalues of the reduced hessian */
    if(eig_rh_eval>0){
        fprintf(stderr, "W[K_AUG]...\t[K_AUG_ASL]""Evaluating the eigenvalues of the solution matrix. \n");
        _is_not_irh_pdf = dsyev_driver(n_dof, x_, K_nrows, hr_point);}
    /*
     * Evaluate the inverse of RedHessian.
     * Note this has been proven unreliable if you have too many bounds active.
     * We need to figure out something better.
     */
    if(comp_inv_ > 0){
        fprintf(stderr, "W[K_AUG]...\t[K_AUG_ASL]""Evaluating the inverse of the solution matrix. \n");
        dpotri_driver(n_dof, x_, K_nrows, hr_point, _chr_timest);}

    total_c = clock();
    /*
     * Quick estimate of our computation time.
     */
    ev_as_time = (double) (ev_as_kkt_c-start_c) / CLOCKS_PER_SEC;
    fact_time = (double) (fact_kkt_c-start_c) / CLOCKS_PER_SEC;
    overall_time = (double) (total_c - start_c) / CLOCKS_PER_SEC;
    printf("I[K_AUG]...\t[K_AUG_ASL]Timings.."
           "Ev&Assem %g, Fact %g, Overall %g.\n",
           ev_as_time, fact_time, overall_time);
    if (dsdp_mode > 0) {
        somefile = fopen("timings_k_aug_dsdp.txt", "a");
    } else {
        somefile = fopen("timings_k_aug.txt", "a");
    }
    /*
     * Get some nice statistics.
     */
    fprintf(somefile, "%d\t%g\t%g\t%g\n", n_var, ev_as_time, fact_time, overall_time);
    fclose(somefile);

    /*
     * Deallocate all pointers.
     */

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
    /*
     * Goodbye.
     */
    return 0;
}
