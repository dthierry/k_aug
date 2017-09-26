#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern void  mc21a_(int *n, int *icn, int *licn, int *ip, int *lenr, int *iperm, int *numnz, int *iw);
extern void mc21ad_(int *n, int *icn, int *licn, int *ip, int *lenr, int *iperm, int *numnz, int *iw);

int main(int argc, char *argv[]){
      int n, m, ne;
      int *icn, *lenr, *ip, *iperm, *iw, licn;
      int *irn, *jcn;
      int i, rn, *rstrt, *rnz;
      int crnz, numnzd;
      float dummy;
      FILE *f_in, *f_out;
      
      if(argc != 2){
            printf("Usage: %s filename_matrix\n", argv[0]);
            return 1; /* terminate if no file is given */
      }
      
      f_in = fopen(argv[1], "r");
      
      fscanf(f_in, "%d\t%d\t%d", &m, &n, &ne);
      printf("%d\t%d\t%d\n", m, n, ne);

      irn = (int *) malloc(sizeof(int)*ne);
      jcn = (int *) malloc(sizeof(int)*ne);
      rstrt = (int *) malloc(sizeof(int)*m);
      rnz   = (int *) malloc(sizeof(int)*m);
      memset(rstrt, 0, sizeof(int) * m);
      memset(rnz, 0, sizeof(int) * m);       

      rn  = 1;
      rstrt[0] = 1;
      crnz = 0;      
      for(i=0; i < ne; i++){
            crnz++;
            fscanf(f_in, "%d\t%d\t%fl", irn + i, jcn + i, &dummy);           
            /*printf("%d\t%d\t%f\n", *(irn + i), *(jcn + i), dummy);*/

            /*foldyflaps*/
            if(rn < *(irn + i)){ /* row changed */
                  rstrt[*(irn + i)-1] = i + 1; /*row_strt*/
                  rn = *(irn + i); /*curent row_number*/
                  rnz[*(irn + i)-2] = crnz; /*update number of non zeroes*/
                  crnz = 0;
            }
            if(*(irn + i) > m){
                  printf("The current row number is somehow greater than the reported number of rows\t");
                  printf("Reported %d\tread %d\n", m, *(irn+i));
                  }

      }

      rnz[*(irn + i-1)-1] = crnz+1;
      fclose(f_in);
      
      f_out = fopen("out_mc21", "w");
      
      for(i=0; i < m; i++){fprintf(f_out, "%d\t%d\n", *(rstrt + i), *(rnz + i));}          
      
      fclose(f_out);
      
      iperm = (int *)malloc(sizeof(int)*m);
      iw = (int *)malloc(sizeof(int)* 5 * m); /*5 just in case */
      mc21ad_(&m, jcn, &ne, rstrt, rnz, iperm, &numnzd, iw);
            
      f_out = fopen("out_permutation.txt", "w");
      for(i=0; i < m; i++) {fprintf(f_out, "%d\n", *(iperm+i));}
      printf("I[[MC21]], number of nz in the permuted diagonal: %d\n", numnzd);
      fclose(f_out);
      
      free(iperm);
      free(iw);
      
      free(irn);
      free(jcn);
      free(rstrt);
      free(rnz);
      return 0;    
}
 
