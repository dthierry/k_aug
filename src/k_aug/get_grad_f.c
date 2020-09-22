
#include "get_grad_f.h"
#include <stdio.h>
#include <stdlib.h>
/*#include "asl.h" */

/* not really necessary, just for consistency */
void get_grad_f(ASL *asl, real *x, fint nvar, real *gf, fint *nerror){
//	gf = (real *)malloc(sizeof(real) * nvar);

	objgrd(0, x, gf, nerror);

	/*for(int i =0; i < nvar; i++){ */
		/*printf("grad_f %f\n", gf[i]); */
	/*} */


//	free(gf);
	/*> tfw no gf */
}
