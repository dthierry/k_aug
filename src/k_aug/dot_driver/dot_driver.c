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

#include "../../../thirdparty/asl/solvers/asl.h"
#include "../../../thirdparty/asl/solvers/getstub.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void dsdp_strategy(ASL *asl, int n_srow, int nvar, SufDesc *dcdp_suf, SufDesc *DeltaP_suf, const char *fname);

void npdp_strategy(ASL *asl, int n_srow, int nvar, SufDesc *suf1, SufDesc *suf2, const char *fname);

extern void
dgemv_(char *TRANS, fint *M, fint *N, real *ALPHA, real *A, fint *LDA, real *X, fint *INCX, real *BETA, real *Y,
       fint *INCY);

static char _dsdp_mode[] = {"dsdp_mode"};
static char _dsdp_verb[] = {"hat{s} = s^{*} + S * DeltaP mode."};

static int dsdp_mode_active = 0;
static I_Known dsdp_mode_active_kw = {1, &dsdp_mode_active};


static keyword keywds[] = {
        KW(_dsdp_mode, IK_val, &dsdp_mode_active_kw, _dsdp_verb)};

static Option_Info Oinfo;


int main(int argc, char **argv) {
    ASL *asl;

    int nvar, ncon;
    int n_srow;

    char *s = NULL;
    FILE *f, *f_out;

    SufDecl *suf_ptr = NULL;

    char _suf1[] = {"dof_v"};
    char _suf2[] = {"npdp"};
    char _suf3[] = {"f_timestamp"};
    char _suf4[] = {"dcdp"};
    char _suf5[] = {"DeltaP"};

    SufDesc *suf1 = NULL;
    SufDesc *suf2 = NULL;
    SufDesc *suf3 = NULL;
    SufDesc *suf4 = NULL;
    SufDesc *suf5 = NULL;


    char ter_msg[] = {"I[[DOT_SENS]]...\t[MAIN]"
                      "All done."};

    char _sname[] = {"[[DOT_SENS]]"};
    clock_t start_c, end_c;
    double cpu_timing;
    time_t timestamp;
    char _chr_timest[15] = "";
    char _file_name_[30] = ""; /*  */


    start_c = clock();
    timestamp = time(NULL);

    Oinfo.sname = _sname;
    Oinfo.bsname = NULL;
    Oinfo.opname = NULL;
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

    asl = ASL_alloc(ASL_read_fg);
    s = getstops(argv, &Oinfo);

    if (!s) {
        printf("I[[DOT_SENS]]...\t[MAIN]"
               "No input.\n");
        ASL_free(&asl);
        return 0;
    } else {
        printf("I[[DOT_SENS]]...\t[MAIN]"
               "File read succesful.\n");
    }

    /* Five suffixes */
    suf_ptr = (SufDecl *) M1alloc(sizeof(SufDecl) * 5);

    /* consider using other suffix instead */
    suf_ptr->name = _suf1;
    suf_ptr->table = 0;
    suf_ptr->kind = ASL_Sufkind_var;
    suf_ptr->nextra = 0;

    (suf_ptr + 1)->name = _suf2;
    (suf_ptr + 1)->table = 0;
    (suf_ptr + 1)->kind = ASL_Sufkind_con | ASL_Sufkind_real;
    (suf_ptr + 1)->nextra = 0;

    (suf_ptr + 2)->name = _suf3;
    (suf_ptr + 2)->table = 0;
    (suf_ptr + 2)->kind = ASL_Sufkind_prob | ASL_Sufkind_input;
    (suf_ptr + 2)->nextra = 0;

    (suf_ptr + 3)->name = _suf4; /* dcdp */
    (suf_ptr + 3)->table = 0;
    (suf_ptr + 3)->kind = ASL_Sufkind_con | ASL_Sufkind_input;
    (suf_ptr + 3)->nextra = 0;

    (suf_ptr + 4)->name = _suf5; /* DeltaP */
    (suf_ptr + 4)->table = 0;
    (suf_ptr + 4)->kind = ASL_Sufkind_con | ASL_Sufkind_input;
    (suf_ptr + 4)->nextra = 0;


    suf_declare(suf_ptr, 5);

    f = jac0dim(s, strlen(s));
    nvar = n_var;
    ncon = n_con;
    n_srow = nvar + ncon;
    printf("I[[DOT_SENS]]...\t[MAIN]"
           "Number of variables %d\n", (int) nvar);
    printf("I[[DOT_SENS]]...\t[MAIN]"
           "Number of constraints %d\n", (int) ncon);
    printf("I[[DOT_SENS]]...\t[MAIN]"
           "Number of rows (primal-dual) %d\n", (int) n_srow);

    /* Primal and dual */
    X0 = M1alloc(sizeof(double) * nvar);
    pi0 = M1alloc(sizeof(double) * ncon);

    fg_read(f, 0);

    suf1 = suf_get(_suf1, ASL_Sufkind_var); /* rh_name */
    suf2 = suf_get(_suf2, ASL_Sufkind_con); /* npdp */
    suf3 = suf_get(_suf3, ASL_Sufkind_prob);
    suf4 = suf_get(_suf4, ASL_Sufkind_prob); /* dcdp starts at 1 instead of 0th */
    suf5 = suf_get(_suf5, ASL_Sufkind_prob); /* DeltaP */
    if (dsdp_mode_active) {
        printf("\n\nI[[DOT_SENS]]...\t[MAIN]"
               "dsdp_mode_active\n\n");
        /*  */
        if (!(suf4->u.i)) {
            printf("E[[DOT_SENS]]...\t[MAIN]"
                   "No \"%s\" suffix declared. Exiting.\n\n", _suf4);
            ASL_free(&asl);
            exit(-1);
        }
        if (!(suf5->u.i)) {
            printf("E[[DOT_SENS]]...\t[MAIN]"
                   "No \"%s\" suffix declared. Exiting.\n\n", _suf5);
            ASL_free(&asl);
            exit(-1);
        }
    } else {
        if (!(suf1->u.i)) {
            printf("E[[DOT_SENS]]...\t[MAIN]"
                   "No \"%s\" suffix declared. Exiting.\n\n", _suf1);
            ASL_free(&asl);
            exit(-1);
        }
        if (!(suf2->u.r)) {
            printf("E[[DOT_SENS]]...\t[MAIN]"
                   "No \"%s\" suffix declared. Exiting.\n\n", _suf2);
            ASL_free(&asl);
            exit(-1);
        }
    }

    if (dsdp_mode_active) {
        strcat(_file_name_, "dsdp_in_");
    } else {
        strcat(_file_name_, "dot_in_");
    }

    if (!(suf3->u.i)) {
        fprintf(stderr, "E[[DOT_SENS]]...\t[MAIN]"
                        "No \"%s\" suffix declared. Fallback to default filename.\n\n", _suf3);
        /*exit(-1);*/
    } else {
        printf("I[[DOT_SENS]]...\t[MAIN]"
               "Timestamp suffix = %d.\n\n", *(suf3->u.i));
        /*_chr_timest[0] = '\0';*/
        sprintf(_chr_timest, "%d", *(suf3->u.i));
        fprintf(stderr, "This goes here %s\n", _chr_timest);
    }
    strcat(_file_name_, _chr_timest);
    strcat(_file_name_, ".in");
    fprintf(stderr, "I[[DOT_SENS]]...\t[MAIN] %s\n", _file_name_);

    /* we branch on strategies
    */
    if (dsdp_mode_active) {
        dsdp_strategy(asl, n_srow, nvar, suf4, suf5, _file_name_);
    } else {
        npdp_strategy(asl, n_srow, nvar, suf1, suf2, _file_name_);
    }

    suf_iput(_suf3, ASL_Sufkind_prob | ASL_Sufkind_iodcl, (int *) &timestamp);
    solve_result_num = 0;
    write_sol(ter_msg, X0, pi0, 0);

    end_c = clock();

    cpu_timing = (double) (end_c - start_c) / CLOCKS_PER_SEC;

    printf("I[[DOT_SENS]]...\t[MAIN]"
           "дава́й!.\tDone.\n");

    printf("I[[DOT_SENS]]...\t[MAIN]Timing.."
           "%g sec.\n", cpu_timing);
    f_out = fopen("timings_dot_driver.txt", "w");
    fprintf(f_out, "%g\n", cpu_timing);
    fclose(f_out);
    ASL_free(&asl);
    return 0;

}

void npdp_strategy(ASL *asl, int n_srow, int nvar, SufDesc *suf1, SufDesc *suf2, const char *fname) {
    int i, j, _n_dof, ncon = (n_srow - nvar), *u_arr = NULL;
    double *npdp = NULL, *u_star = NULL, *u_star0 = NULL;
    FILE *rh_txt = NULL, *f_out = NULL;
    double *s_hat_T = NULL;

    char t = 'T';
    double ALPHA = -1.0;
    /*int LDA = 10;*/
    int INCX = 1;
    double BETA = 1.0;
    int INCY = 1;

    npdp = (double *) malloc(sizeof(double) * n_srow);
    memset(npdp, 0, sizeof(double) * n_srow);

    for (i = 0; i < ncon; i++) {
        *(npdp + nvar + i) = suf2->u.r[i];
    }

    u_arr = (int *)malloc(sizeof(int) * nvar); /* This guy points to the E^T s*. To the required rows to be precise */
    memset(u_arr, 0, sizeof(int) * nvar);
    _n_dof = 0; /* count_dof */
    for (i = 0; i < nvar; i++) {
        if (*((suf1->u.i) + i) != 0) {
            u_arr[_n_dof] = (int) i;
            _n_dof++;
        }
    }

    u_star = (double *) calloc(_n_dof, sizeof(double)); /* The primal solution */
    u_star0 = (double *) malloc(sizeof(double) * _n_dof);
    /* Load the E^T s* primal*/
    for (i = 0; i < _n_dof; i++) {
        u_star[i] = asl->i.X0_[u_arr[i]];
        u_star0[i] = u_star[i];
    }
    printf("I[[DOT_SENS]]...\t[NPDP_STRATEGY]"
           "Number of dof detected  %u\n", _n_dof);
    rh_txt = fopen(fname, "r");
    /* Flat file with solution vectors of Shat*/
    if (!rh_txt) {
        fprintf(stderr, "W[[DOT_SENS]]...\t[MAIN]File %s not found.\n", fname);
        ASL_free(&asl);
        free(u_arr);
        free(u_star);
        exit(-1);
    }

    s_hat_T = (double *) malloc(sizeof(double) * n_srow * _n_dof);

    for (i = 0; i < _n_dof; i++) {
        for (j = 0; j < n_srow; j++) {
            fscanf(rh_txt, "%lf", (s_hat_T + n_srow * i + j));
        }
    }
    fclose(rh_txt);
    dgemv_(&t, (fint *) &n_srow, (fint *) &_n_dof, &ALPHA, s_hat_T, (fint *) &n_srow, npdp, &INCX, &BETA, u_star,
           &INCY);
    /*sens_update_driver_dot((fint) n_srow, (fint) _n_dof, s_hat_T, npdp, u_star);*/

    f_out = fopen("dot_out.out", "w");

    for (i = 0; i < _n_dof; i++) {
        fprintf(f_out, "%.g\t%.g\t%.g\n", u_star0[i], u_star[i], u_star0[i] - u_star[i]);
    }
    fclose(f_out);


    for (i = 0; i < _n_dof; i++) { asl->i.X0_[u_arr[i]] = u_star[i]; } /* update primal */

    free(u_star);
    free(u_star0);
    free(u_arr);
    free(npdp);
    free(s_hat_T);
}

void dsdp_strategy(ASL *asl, int n_srow, int nvar, SufDesc *dcdp_suf, SufDesc *DeltaP_suf, const char *fname) {

    int i, j, _n_p = 0, ncon = (n_srow - nvar), temp;
    double *dpvect = NULL, *s_star0 = NULL, *s_star = NULL;
    FILE *s_txt = NULL;
    double *s_ = NULL;
    char t = 'N';
    double ALPHA = -1.0;
    /*int LDA = 10;*/
    int INCX = 1;
    double BETA = 1.0;
    int INCY = 1;


    dpvect = (double *) calloc(n_srow, sizeof(double));
    s_star = (double *) calloc(n_srow, sizeof(double));
    s_star0 = (double *) calloc(n_srow, sizeof(double));

    /* load solution */
    for (i = 0; i < nvar; i++) {
        s_star[i] = asl->i.X0_[i];
        s_star0[i] = s_star[i];
    }

    _n_p = 0;
    /* load deltaP*/
    /* dcdp contains; i.e. the order */
    for (i = 0; i < ncon; i++) {
        temp = dcdp_suf->u.i[i]; /* Retrieve value of suffix*/
        /* Find non-zero */
        if (temp != 0) {
            if (temp - 1 >= ncon) {/* error */
                printf("E[DOT_SENS]...\t[.]"
                       "The suffix dcdp at %d is greater than n_con e.g. %d\n", i, ncon);
                exit(-1);

            } else if (temp - 1 < 0) {/* error again*/
                printf("E[DOT_SENS]...\t[.]"
                       "The suffix dcdp at %d is negative (%d)\n", i, temp);
                exit(-1);

            }
            /* put to dpvect temp starts at 1*/
            dpvect[temp - 1] = DeltaP_suf->u.r[i];
            _n_p++;
        }
    }

    /* */
    printf("I[[DOT_SENS]]...\t[DSDP_STRATEGY]"
           "Number of parameters detected  %d\n", _n_p);

    s_txt = fopen(fname, "r");
    /* Flat file with solution vectors of s*/
    if (!s_txt) {
        fprintf(stderr, "W[[DOT_SENS]]...\t[DSDP_STRATEGY]File %s not found.\n", fname);
        ASL_free(&asl);
        free(dpvect);
        free(s_star0);
        free(s_star);
        exit(-1);
    }
    /* sensitivity matrix */
    s_ = (double *) malloc(sizeof(double) * n_srow * _n_p);

    for (i = 0; i < _n_p; i++) {
        for (j = 0; j < n_srow; j++) {
            fscanf(s_txt, "%lf", (s_ + n_srow * i + j));
        }
    }
    fclose(s_txt);

    dgemv_(&t, (fint *) &n_srow, (fint *) &_n_p, &ALPHA, s_, (fint *) &n_srow, dpvect, &INCX, &BETA, s_star, &INCY);
    /* */
    for (i = 0; i < n_var; i++) { asl->i.X0_[i] = s_star[i]; } /* update primal */

    free(dpvect);
    free(s_star0);
    free(s_star);
    free(s_);
}