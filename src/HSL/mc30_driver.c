/* @source mc30_driver.c
**
** April 10th, 2017
** @author: David Thierry (dmolinat@andrew.cmu) dav0@lb2016-1

******************************************************************************

@mc30driver ********************************************
**
** Calculates scaling factors
**
** @param [r] nr row number
** @param [r] nz number of nz
** @@
*******************************************************************************/

/*#include <stdlib.h> */
#include "mc30_driver.h"

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>


int mc30driver(fint n, fint nz, real *a, fint *irn, fint *icn, real *s){
	
	real *w;
	fint ifail, lp;
	FILE *somefile;
	int i;
	w = (real *)malloc(sizeof(real) * n * 4);
	memset(w, 0.0, n*4);
	memset(s, 0.0, n);
	
	/* initialize s & w */
	/*for(i=0; i<n; i++){ */
	/*	s[i] = 0.0; */
	/*} */

	/*for(i=0; i<(4*n); i++){ */
	/*	w[i] = 0.0; */
	/*} */

	ifail = 0;
	lp = 6;

	/* call mc30ad */
	mc30ad_(&n, &nz, a, irn, icn, s, w, &lp, &ifail);
	
	assert(ifail == 0);
	free(w);

	somefile = fopen("scale_facts.txt", "w");
	for(i=0; i < n; i++){
		fprintf(somefile, "%f\n", s[i]);
	}
	fclose(somefile);

	return 0;
}