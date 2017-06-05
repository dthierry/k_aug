/* @source untitled.c
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

#include "asl.h"
#include "getstub.h"
#include <stdio.h>
#include <stdlib.h>
static Option_Info Oinfo;
static keyword *keywds;

int main(int argc, char **argv){
	ASL *asl;
	unsigned i, j, k;
	char *s=NULL;
	FILE *f;

	double *x, *lambda;

  char ter_msg[] = {"I[[SENS_DOT]]...\t[MAIN]"
	"All done it seems."};

	char _sname[] = {"[[SENS_DOT]]"};

	Oinfo.sname = _sname;
	Oinfo.bsname = NULL;
	Oinfo.opname = NULL;
	Oinfo.keywds = NULL;
	Oinfo.n_keywds = nkeywds;
	Oinfo.flags = 0;
	Oinfo.version = NULL;
	Oinfo.usage = NULL;
	Oinfo.kwf = NULL;
	Oinfo.feq = NULL;
	Oinfo.n_options = 0;
	Oinfo.driver_date = 0;
	Oinfo.wantsol = 0;
	Oinfo.nS = 0;
	Oinfo.S = NULL;
	Oinfo.uinfo = NULL;
	Oinfo.asl = NULL;
	Oinfo.eqsign = NULL;
	Oinfo.n_badopts = 0;
	Oinfo.option_echo = 0;
	Oinfo.nnl = 0;

	asl = ASL_alloc(ASL_read_fg);
	s = getstops(argv, &Oinfo);

	if(!s){
		printf("I[[SENS_DOT]]...\t[MAIN]"
			"No input.\n");
		ASL_free(&asl);
		return 0;
	}
	else{
		printf("I[[SENS_DOT]]...\t[MAIN]"
			"File read succesful.\n");
	}

	f = jac0dim(s, (int)strlen(s));
	printf("nvar %d ncon %d\n", n_var, n_con);

	/* Primal and dual */
	X0 = M1alloc(sizeof(double) * n_var);
	pi0= M1alloc(sizeof(double) * n_con);

	fg_read(f, 0);

	x = X0;
	lambda = pi0;

	ASL_free(&asl);
	return 0;

}