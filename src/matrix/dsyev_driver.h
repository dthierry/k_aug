/* @source untitled.c
** beta 01
** Month dayth, 20yy
** @author: David Thierry (dmolinat@andrew.cmu) dav0@lb2016-1

********************************************************************************

@fun_name ********************************************
**
** Driver for eigenvalue calculation
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


#ifndef DSYEV_DRI
#define DSYEV_DRI
#include "../../thirdparty/asl/solvers/asl.h"
extern void dsyev_(char *jobz, char* uplo, long *n, double *a,
	int *lda, double *w, double *work, int *lwork, int *info);


int dsyev_driver(long n, double *a, long Kn, int *sb_p);

#endif
