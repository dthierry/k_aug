#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern void mc21a_(int *n, int *icn, int *licn, int *ip, int *lenr, int *iperm, int *numnz, int *iw);
extern void mc21ad_(int *n, int *licn, int *ip, int *iperm, int *numnz, int *iw);

int main(int argc, char *argv[]){
      int n, m, ne;
      int *icn, *lenr, *ip, *iperm, *iw, licn;
      int *irn, *jcn;
      int i, rn, *rstrt, *rnz;
      int crnz;
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
                  rnz[*(irn + i)-1] = crnz; /*update number of non zeroes*/
                  crnz = 0;
            }
            if(*(irn + i) > m){
                  printf("The current row number is somehow greater than the reported number of rows\t");
                  printf("Reported %d\tread %d\n");
                  }
      } 
      
      fclose(f_in);
      
      f_out = fopen("out_mc21", "w");
      
      for(i=0; i < m; i++){fprintf(f_out, "%d\t%d\n", *(rstrt + i), *(rnz + i));}          
      
      fclose(f_out);
      free(irn);
      free(jcn);
      free(rstrt);
      free(rnz);
      return 0;    
}
 
