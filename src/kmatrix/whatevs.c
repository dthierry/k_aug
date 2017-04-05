#include "getstub.h"
 static int dumm = 1;
 static I_Known dumm_kw = {2, &dumm};
static int n_rhs = 0;


  // suffix
  static SufDecl
  suftab[] = {
      { "var_flag", 0, ASL_Sufkind_var, 0},
      { "con_flag", 0, ASL_Sufkind_con, 0},
  };

SufDecl *Suftabdecl(int n_rhs){
  //char rhs[] = "rhs_";
  //int size = 2;
  SufDecl *s;
  if(n_rhs == 0){
    printf("Warning: no rhs_n declared");
  }
  else{
    printf("Number of rhs for suftab %d\n", n_rhs);
  }
  SufDecl suftabx[n_rhs + 2];
  s = suftabx;
  int k = sizeof(SufDecl);
  printf("The size of sufdecl is %d\n", k);
  return s;
};

 // keywords
static keyword keywds[] = {
  KW("smth", IK_val, &dumm_kw, "Cheers mate the cavalry is here"),
  KW("n_rhs", I_val, &n_rhs, "Number or right hand sides")
  };

 static Option_Info Oinfo = 
{"whatevs", "whatevs", "whatevs_options", keywds, nkeywds};
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
      SufDesc *var_f;
      SufDesc *con_f;
      SufDecl *tabx;
//      real ow;
      /*Pointer to file name*/
      int *hcolstrt, *hrown;
  
	asl = ASL_alloc(ASL_read_pfgh);
	/* The memory allocation */

      s = getstops(argv, &Oinfo);

      if (!s) {
            return 1;
            }
      else {
            printf("File read succesful\n");
            }
      tabx = Suftabdecl(n_rhs);

      suf_declare(suftab, sizeof(suftab)/sizeof(SufDecl));

      f = jac0dim(s,0);
      
      printf("Number of variables \t%d\n", n_var);
      printf("Number of constraints \t%d\n", n_con);     

      x = X0 = M1alloc(n_var*sizeof(real));
      
      // Want initial guess
      want_xpi0 = 2;
      
      // need to do part of changing sign for y 
      
      pfgh_read(f, ASL_findgroups);
 
      y = pi0;

      var_f = suf_get("var_flag", ASL_Sufkind_var);
      con_f = suf_get("con_flag", ASL_Sufkind_con);

      //printf("Whatch out!!%f\n", zLmult->u.r);

      //zLmult = zLmult->next;

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
 
     // printf("Whatch out!!%d\n", zLmult->u.r);
     //      for(i = 0; i < n_con; i++){
     //            printf("pointer array at %d y %f\n", i, *(y + i));
     //            };

      nh = sphsetup(0, 0, 1, 0);
      // setup the sparse hessian with multipliers
      i = strlen(filename) + 60;

      printf("i strlen() + 60 %d\n", i);
      for(i=0; i < strlen(filename); i++){
            printf("%c\n",filename[i]);
      }

      f_jac = fopen("jacobi_debug.in", "w");
      // Little file for the dummy results    

      // Assert allocation factor for gradient of the objective
      nn = 100;
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


      return 0;

}
