/* @source dpotri_driver.c
** beta 01
** June 26th, 2017
** @author: David Thierry (dmolinat@andrew.cmu) dav0@lb2016-1

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
#ifndef DPOTRI_DRIVER
#define DPOTRI_DRIVER

extern void dpotrf_(char *uplo, int *n, double *a, int *lda, int *info);
extern void dpotri_(char *uplo, int *n, double *a, int *lda, int *info);

void dpotri_driver(int n, double *_a, long Kn, int *sb_p, char *_chr_timest);

#endif
