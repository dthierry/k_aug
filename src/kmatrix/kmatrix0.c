#include "getstub.h"
#include <assert.h>
#include "csortv2.h"
#include "wapnz.h"
#include "wcreord.h"


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
        FILE *f_debug;
        /* SufDesc *some_suffix; */
        int i, j, k;
        int nnzw, nn; // let's try this
        cgrad *cg, **cgx;
        real *x, *y;
        real *g, *g1;
        char *s;
        static char sword[] = "rhs_";
        SufDesc *var_f = NULL;
        SufDesc *con_f = NULL;
        SufDesc *rhs_ptr = NULL;
        SufDecl *tptr;
        // Pointer to file name
        fint *Acol, *Arow;
        real *Aij;
        fint *Wcol, *Wrow;
        real *Wij;

        // The objective weight
        real ow;

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

        // (tptr + 2)->name = "rhs_0";
        // (tptr + 2)->table = 0; 
        // (tptr + 2)->kind = ASL_Sufkind_con;
        // (tptr + 2)->nextra = 0;

        // (tptr + 3)->name = "rhs_1";
        // (tptr + 3)->table = 0; 
        // (tptr + 3)->kind = ASL_Sufkind_con;
        // (tptr + 3)->nextra = 0;

        // Suffix for the right hand side creation
        
        //char *suffix_name[] = {"rhs_0", "rhs_1"};
        char **suf_name;
        suf_name = (char **)malloc(sizeof(char *)*n_rhs);
       
        for(i=0; i < n_rhs; i++){
                char numeric_rhs[5];
                numeric_rhs[0] = '\0'; 
                //*(suf_name[i]) = '\0';
                sprintf(numeric_rhs, "%d", i);
                suf_name[i] = (char *)malloc(sizeof(char)*10);
                *(suf_name[i]) = '\0';

                strcat(suf_name[i], sword);
                strcat(suf_name[i], numeric_rhs);
                printf("The name of the suffix %s\n", suf_name[i]); 
        }       
        
        for(i=0; i<n_rhs; i++){
        printf("NAME %s\n", suf_name[i]);
        }

        for(i=0; i<n_rhs; i++){
                (tptr + i + 2)->name = suf_name[i];
                (tptr + i + 2)->table = 0; 
                (tptr + i + 2)->kind = ASL_Sufkind_con;
                (tptr + i + 2)->nextra = 0;
                //free(suf_name);
        }

        //rhs_ptr = (SufDesc **)malloc(sizeof(SufDesc *)*2);

        // Declare suffixes
        suf_declare(tptr, (n_rhs + 2));
        

        f = jac0dim(s, 0);

        x = X0 = M1alloc(n_var*sizeof(real));

        // Want initial guess
        want_xpi0 = 2;

        // need to do part of changing sign for y 

        pfgh_read(f, ASL_findgroups);

        y = pi0;
        
        con_f = suf_get("con_flag", ASL_Sufkind_con);
        var_f = suf_get("var_flag", ASL_Sufkind_var);

        if(var_f->u.i == NULL){
                fprintf(stderr, "u.i empty, no var_flag declared.\n");
                return -1;
        }
        if(con_f->u.i == NULL){
                fprintf(stderr, "u.i empty, no con_flag declared.\n");
                return -1;
        }

        //SufDesc *smptr0, *smptr1;
        
        // smptr0 = suf_get("rhs_0", ASL_Sufkind_con);
        // smptr1 = suf_get("rhs_1", ASL_Sufkind_con);

        for(i=0; i < n_rhs; i++){
                printf("The suffix name %s \n", suf_name[i]);
                rhs_ptr = suf_get(suf_name[i], ASL_Sufkind_con);
                //if(rhs_ptr[i]->u.r == NULL){
                        //fprintf(stderr, "u.r empty, no rhs values declared for rhs_%d.\n", i);
                //return -1;
                }
        // }


        
        //printf("Pointer address current \t%p\n", var_f);
        //printf("Pointer address next \t%p\n", var_f->next);

        //for(i=0; i < n_con; i++){
                //if(con_f->u.i[i] != 0){
                        //printf("suffix value for i = %d, %d\n",i, con_f->u.i[i]);
                //}
        //}
        //printf("Pointer address current \t%p\n", con_f);
        //printf("Pointer address next \t%p\n", con_f->next);

        //};
        

        // setup the sparse hessian with multipliers
        if (n_obj == 0){
                ow = 0; // set the objective weight to zero
                nnzw = sphsetup(0, ow, 1, 1);
                printf("No objective declared \n");
        }
        else{
                ow = 1; // set the objective weight to one
                nnzw = sphsetup(-1, ow, 1, 1);
                printf("Objective found set ow = 1\n");
        }

        if (n_obj > 0){
                if(objtype[0]){
                        printf("Maximization problem detected\n");
                        // set weight to -1
                        ow = -1;
                }
                else{
                        printf("Minimization problem detected\n");
                        ow = 1;
                }
                //objtype is an int?
                //printf("Number of objectives %d\n", n_obj);
                //printf("%d\n", objtype[0]);
        }

        f_jac = fopen("jacobi_debug.in", "w");
        // Little file for the dummy results

        // Assert allocation factor for gradient of the objective
        nn = 2*n_var;
        if (nn < nnzw)
                nn = nnzw;
        if (nn < nzc)
                nn = nzc;

        g = (real *)M1alloc(nn*sizeof(real));
        g1 = 0;
        if (n_obj > 0){
                // Unless somebody does something wrong there should be only
                // one objective
                objgrd(0, x, g1 = g, 0);
        }
        
        // Row and colum for the triplet format A matrix
        // size of the number of nonzeroes in the constraint jacobian
        Acol = (fint *)malloc(sizeof(fint)*nzc);
        Arow = (fint *)malloc(sizeof(fint)*nzc);
        Aij = (real *)malloc(sizeof(real)*nzc);
        j = 0;
        // Jacobian
        // Start row in the triplet matrix at 0, Actual matrix is in Fortran format
        if (nzc > 0) {
                jacval(x, g, 0);
                cgx = Cgrad; // pointer magic; assuming there is no zero column
                for(i = 1; i <= n_con; i++) {
                        // moves by constraint
                        if ((cg = *cgx++)) {
                                // iterates for a given constraint
                                do{
                                        // moves by nz in the constraint
                                        fprintf(f_jac, "%d\t%d\t%.g\n",i , cg->varno+1, g[cg->goff]);
                                        Arow[j] = i;
                                        Acol[j] = cg->varno+1;
                                        Aij[j] = g[cg->goff];
                                        j++;
                                        }
                                while ((cg = cg->next)) ;
                        }
                }
        };

        fclose(f_jac);

        f_hess = fopen("hess_debug.in", "w");

        
        Wij = (real *)malloc(sizeof(real)*nnzw);
        Wcol = (fint *)malloc(sizeof(fint)*nnzw);
        Wrow = (fint *)malloc(sizeof(fint)*nnzw);
        // Hessian of the Lagrange function matrix
        if (nnzw) {
                if (n_obj > 0){
                        sphes(g, -1, &ow, y);
                }
                else{
                        sphes(g, 0, &ow, y);
                }
                // pretty much compressed column format
                k = 0;
                // position the counter at 0
                for(i = 0; i < n_var; i++){
                        for (j = sputinfo->hcolstarts[i]; j< sputinfo->hcolstarts[i+1]; j++){
                                Wij[j] = g[j];
                                Wcol[j] = i + 1;
                                Wrow[j] = sputinfo->hrownos[j] + 1;
                                fprintf(f_hess, "\t%ld\t%ld\t%.g\n", sputinfo->hrownos[j] + 1, i+1, g[j]);
                                k++;
                        }
                }

        }
        fclose(f_hess);

        printf("nonzeroes in the sparse hessian %d\n", nnzw);
        printf("print dummy %d\n", dumm);
        if (dumm == 2){
                printf("Cheers mate, the cavalry is here!");
        }
        printf("size of s %d\n", sizeof(f));
        printf("Number of Right hand sides %d\n", n_rhs);
        printf("Number of variables %d\n", n_var);
        printf("Number of constraints %d\n", n_con);
        fint nzA = nzc;
        /*
        sortcol(Arow, Acol, Aij, nzA, n_var);
        f_debug = fopen("somefile.txt", "w");
        for(i=0; i < nzA; i++){
                fprintf(f_debug, "%d\t%d\t%.g\n", Arow[i], Acol[i], Aij[i]);
        }
        fclose(f_debug);
        */

        wnzappnd(Wrow, Wcol, Wij, nnzw, n_var);
        //w_col_sort(Wrow, Wcol, Wij, nnzw, n_var);

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
        //free(rhs_ptr);
        return 0;

};
;