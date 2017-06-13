CC = gcc
#CCFLAGS = -g -Wall -ansi -pedantic-errors -Wall -W -Wconversion -Wshadow -Wcast-qual -Wwrite-strings
CCFLAGS = -g -Wall -ansi -pedantic-errors -Wall -W
FC = gfortran
S = ../..
L = -lm -ldl
A = $S/ASL/solvers/amplsolver.a
PARENT = ..
ONAME = k_


all1 = k_aug
MATRIXDIR = ../matrix
LMATRIX = kassemble.f

# this line !!!
MC30ROOT = /home/dav0/Apps/mc30-1.0.0/build/src/
# this line !!!
PARDISOROOT = /usr/local/lib/


#OPENBLASROOT = /home/dav0/Apps/OpenBLAS/


BLASROOT = $(S)/BLAS/BLAS-3.7.0/

LAPACKROOT = $(S)/LAPACK/lapack-3.7.0/


#LMKLIB= -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl   -m64 -I${MKLROOT}/include
LIBPARDISO = -L${PARDISOROOT} -lpardiso500-GNU481-X86-64 -lgfortran -fopenmp -lpthread -lm
# example LIBOPENBLAS= -L${OPENBLASROOT} -lopenblas_haswellp-r0.2.20.dev
#LIBOPENBLAS= -L${OPENBLASROOT} -lopenblas_haswellp-r0.2.20.dev
LIBBLAS = $(BLASROOT)blas_LINUX.a
LIBLAPACK = $(LAPACKROOT)liblapack.a

LIBMC30 = ${MC30ROOT}libmc30.a -lgfortran -lm -lquadmath
KMATRIX = iteration4_main1.c


k_aug_blas: mc30_dri \
      get_jac_asl get_hess_asl find_ineq get_grad_f assemble_rhsds_rh \
      suffix_decl_hand csr_driver pard_dri sig_comp assemble_corrector_rhs
	$(CC) $(CCFLAGS) $(KMATRIX) -o $(ONAME)  $A $L \
	mc30_driver.o  \
	get_jac_asl_aug.o get_hess_asl_aug.o find_inequalities.o get_grad_f.o \
	assemble_rhs_rh.o  \
	suffix_decl_hand.o csr_driver.o pardiso_driver.o \
	sigma_compute.o assemble_corrector_rhs_v2.o \
	-I$(S)/HSL/ -I$(S)/Pardiso/ -I$(S)/matrix/ -I$(PARENT) -I$(S)/ASL/solvers/\
	 $(LIBMC30) $(LIBPARDISO) $(LIBLAPACK) $(LIBBLAS)


clean:
	rm -rf $(all) *.o
      
mc30_dri:
	$(CC) $(CCFLAGS) -c $(S)/HSL/mc30_driver.c -I$(S)/HSL/ -I$S/ASL/solvers/

pard_dri:
	$(CC) $(CCFLAGS) -c $(S)/Pardiso/pardiso_driver.c -I$(S)/ASL/solvers/


get_jac_asl:
	$(CC) $(CCFLAGS) -c get_jac_asl_aug.c -I$S/ASL/solvers/

get_hess_asl:
	$(CC) $(CCFLAGS) -c get_hess_asl_aug.c -I$S/ASL/solvers/

find_ineq:
	$(CC) $(CCFLAGS) -c $(PARENT)/find_inequalities.c -I$S/ASL/solvers/

get_grad_f:
	$(CC) $(CCFLAGS) -c $(PARENT)/get_grad_f.c -I$S/ASL/solvers/

assemble_corrector_rhs:
	$(CC) $(CCFLAGS) -c assemble_corrector_rhs_v2.c -I$S/ASL/solvers/ -I$(PARENT)

assemble_rhsds_v3:
	$(CC) $(CCFLAGS) -c $(PARENT)/assemble_rhsdsv3.c -o assemble_rhsds.o -I$(S)/ASL/solvers/

sens_update_driver:
	$(CC) $(CCFLAGS) -c $(S)/matrix/sens_update_driver.c -I$S/ASL/solvers/

assemble_rhsds_rh:
	$(CC) $(CCFLAGS) -c assemble_rhs_rh.c -I$(S)/ASL/solvers/

suffix_decl_hand:
	$(CC) $(CCFLAGS) -c suffix_decl_hand.c -I$(S)/ASL/solvers/

csr_driver:
	$(CC) $(CCFLAGS) -c csr_driver.c -I$(S)/ASL/solvers/

sig_comp:
	$(CC) $(CCFLAGS) -c sigma_compute.c -I$(S)/ASL/solvers/
