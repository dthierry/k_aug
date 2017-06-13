/* @source assemble_kkt_by_row.
** beta 01
** Month dayth, 20yy
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
#ifndef CSR_DRIVER
#define CSR_DRIVER

#include <asl.h>


void csr_driver(int nvar, int ncon, int nzW, int nzA,
	int *nzr_w, int *nzr_a,
	int *Wrow, int *Wcol, real *Wij,
	int *Arow, int *Acol, real *Aij,
	fint **Kr, fint **Kc, real **K, fint **Kr_strt);

#endif /* CSR_DRIVER */