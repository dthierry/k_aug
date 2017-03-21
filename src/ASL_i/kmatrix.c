#include "getstub.h"
static int dumm = 1;
static I_Known dumm_kw = {2, &dumm};
static int n_rhs = 0;

// keywords
static keyword keywds[] = {
        KW("smth", IK_val, &dumm_kw, "Cheers mate the cavalry is here"),
        KW("n_rhs", I_val, &n_rhs, "Number or right hand sides")
};

static Option_Info Oinfo = 
{"k_matrix", "KKT_matrix_man", "k_matris_opt", keywds, nkeywds};
/*Remember to set reference back*/

int main(int argc, char **argv)
{
        ASL *asl;
        FILE *f, *f_jac, *f_hess;
        /* SufDesc *some_suffix; */
        int i, ii, iii, j;
        int nh, nn; // let's try this
        cgrad *cg, **cgx;
        real *x, *y;
        real *g, *g1;
        char *s;
        static char sword[] = "rhs_";
        SufDesc *var_f;
        SufDesc *con_f;
        SufDecl *tptr;
        // Pointer to file name
        int *hcolstrt, *hrown;
        
        // The memory allocation
        asl = ASL_alloc(ASL_read_pfgh);

        
        s = getstops(argv, &Oinfo);

        if (!s) {
                return 1;
        }
        else {
                printf("File read succesful\n");
        }

        if (n_rhs == 0){
                fprintf(stderr, "No right hand sides declared\n");
                return -1;
        }


        printf("Number of Right hand sides %d\n", n_rhs);

        // Declare suffix array pointer
        tptr = (SufDecl *)malloc(sizeof(SufDecl)*(2 + n_rhs));

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
        for(i=0; i < n_rhs; i++){
                char numeric_rhs[10];
                char suf_name[16];
                suf_name[0] = '\0';
                sprintf(numeric_rhs, "%d", i);
                printf("Conversion to string %s\n", numeric_rhs);
                strcat(suf_name, sword);
                strcat(suf_name, numeric_rhs);
                printf("Suffix name string %s\n", suf_name);
                (tptr + i + 2)->name = suf_name;
                (tptr + i + 2)->table = 0; 
                (tptr + i + 2)->kind = ASL_Sufkind_con;
                (tptr + i + 2)->nextra = 0;
        }

        // Declare suffixes
        suf_declare(tptr, 2 + n_rhs);

        f = jac0dim(s, 0);

        x = X0 = M1alloc(n_var*sizeof(real));

        // Want initial guess
        want_xpi0 = 2;

        // need to do part of changing sign for y 

        pfgh_read(f, ASL_findgroups);

        y = pi0;

        var_f = suf_get("var_flag", ASL_Sufkind_var);
        con_f = suf_get("con_flag", ASL_Sufkind_con);



        for(i=0; i < n_var; i++){
                if(var_f->u.i[i] != 0){
                        printf("suffix value for i = %d, %d\n",i, var_f->u.i[i]);
                }
        }
        printf("Pointer address current \t%p\n", var_f);
        printf("Pointer address next \t%p\n", var_f->next);

        for(i=0; i < n_con; i++){
        if(con_f->u.i[i] != 0){
        printf("suffix value for i = %d, %d\n",i, con_f->u.i[i]);
        }
        }
        printf("Pointer address current \t%p\n", con_f);
        printf("Pointer address next \t%p\n", con_f->next);


        nh = sphsetup(0, 0, 1, 0);
        // setup the sparse hessian with multipliers

        f_jac = fopen("jacobi_debug.in", "w");
        // Little file for the dummy results    

        // Assert allocation factor for gradient of the objective
        nn = 2*n_var;
        if (nn < nh)
                nn = nh;
        if (nn < nzc)
                nn = nzc;

        g = (real *)M1alloc(nn*sizeof(real));
        g1 = 0;
        if (n_obj > 0){
                // Unless somebody does something wrong there should be only
                // one objective
                objgrd(0, x, g1 = g, 0);
        }
        //    Skip the rest of the useless stuff that I don't need
        //    Column wise accordingly
        if (nzc > 0) {
                jacval(x, g, 0);
                cgx = Cgrad; // pointer magic; assuming there is no zero colmn
                for(i = 1; i <= n_con; i++) {
                        if ((cg = *cgx++)) {
                                do 
                                        fprintf(f_jac, "%d\t%d\t%.g\n",i , cg->varno+1, g[cg->goff]);
                                while ((cg = cg->next)) ;
                        }
                }
        };
        fclose(f_jac);

        f_hess = fopen("hess_debug.in", "w");

        if (nh) {
                sphes(g, 0, NULL, y);
                hcolstrt = sputinfo->hcolstarts;
                hrown = sputinfo->hrownos;
                // pretty much compressed column format

                // position the counter at 0
                ii = 0;
                iii = 0;      
                for(j=1; j<=n_var; j++){
                        iii = hcolstrt[j];
                        if (ii < iii) {
                                do
                                        fprintf(f_hess, "\t%ld\t%ld\t%.g\n", hrown[ii] + 1, g[ii]);
                                while(++ii < iii);
                        }
                }
        }
        fclose(f_hess);

        printf("nonzeroes in the sparse hessian %d\n", nh);
        printf("print dummy %d\n", dumm);
        if (dumm == 2){
                printf("Cheers mate, the cavalry is here!");
        }
        printf("size of s %d\n", sizeof(f));
        printf("Number of Right hand sides %d\n", n_rhs);
        ASL_free(&asl);
        free(tptr);

        return 0;

}
