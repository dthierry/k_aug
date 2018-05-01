/*ExampleprogramusingtheCinterfacetothe∗doubleprecisionversionofMUMPS,dmumpsc.∗WesolvethesystemAx=RHSwith∗A=diag(12)andRHS=[14]ˆT∗Solutionis[12]ˆT*/

#include <stdio.h>
#include "mpi.h"
#include "dmumps_c.h"

#define JOBINIT -1
#define JOBEND -2
#define USECOMMWORLD -987654

intmain(int argc, char **argv){

DMUMPS_STRUC_C id;
int n=2;
int64_t nnz=2;
int jcn[]={1,2};
int irn[]={1,2};
double a[2];
double rhs[2];
int myid, ierr;

ierr = MPIInit(&argc, &argv);
ierr = MPICommrank(MPICOMMWORLD,&myid);

/*DefineAandrhs*/

rhs[0]=1.0;
rhs[1]=4.0;
a[0]=1.0;
a[1]=2.0;

/*InitializeaMUMPSinstance.UseMPICOMMWORLD.*/

id.job=JOBINIT;
id.par=1;
id.sym=0;
id.commfortran=USECOMMWORLD;

dmumpsc(&id);
/*Definetheproblemonthehost*/
if(myid==0){
	id.n=n;
	id.nnz=nnz;
	id.irn=irn;
	id.jcn=jcn;
	id.a=a;
	id.rhs=rhs;
}
#define ICNTL(I) icntl[(I)-1] /*macros.t.indicesmatchdocumentation*/
/*Nooutputs*/

id.ICNTL(1)=-1;
id.ICNTL(2)=-1;
id.ICNTL(3)=-1;
id.ICNTL(4)=0;
/*CalltheMUMPSpackage.*/

id.job=6;
dmumps_c(&id);

id.job=JOBEND;

dmumps_c(&id);/*Terminateinstance*/

if(myid==0){
	printf("Solution is:(%8.2f%8.2f)\n",rhs[0],rhs[1]);
}

ierr = MPI_Finalize();

return 0;
}
