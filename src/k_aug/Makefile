CC = gcc
#CCFLAGS = -g -Wall -ansi -pedantic-errors -Wall -W -Wconversion -Wshadow -Wcast-qual -Wwrite-strings
CCFLAGS = -g -Wall -ansi -pedantic-errors -Wall -W
FC = gfortran
S = ..
ROOT = ../../
L = -lm -ldl
A = $S/ASL/solvers/amplsolver.a
ONAME = k_aug


all1 = k_aug
MATRIXDIR = ../matrix
LMATRIX = kassemble.f

mkfile_path := $(abspath $(lastword $(MAKEFILE_LIST)))
current_dir := $(patsubst %/,%,$(dir $(mkfile_path)))

OPENBLASLIB = -L$(ROOT)thirdparty/openblas/OpenBLAS/ -lopenblas
MUMPSLIB = -L$(current_dir)/../../thirdparty/mumps/MUMPS/lib/ -ldmumps -lmumps_common
MUMPSINC = -I/$(current_dir)/../../thirdparty/mumps/MUMPS/include/
ASLLIB = $(ROOT)src/ASL/solvers/amplsolver.a
ASLINC = -I$(current_dir)/../src/ASL/solvers/
LIBSEQ = -L$(current_dir)/../../thirdparty/mumps/MUMPS/libseq/ -lmpiseq
INCSEQ = -I$(current_dir)/../../thirdparty/mumps/MUMPS/libseq/

LSCOTCHDIR = $(current_dir)/../../thirdparty/scotch/scotch/lib/
LSCOTCH   = -L$(LSCOTCHDIR) -lesmumps -lscotch -lscotcherr
LMETISDIR = $(current_dir)/../../thirdparty/metis/metis/build/Linux-x86_64/libmetis/
LMETIS    = -L$(LMETISDIR) -lmetis

LPORDDIR = -L$(current_dir)/../../thirdparty/mumps/MUMPS/PORD/lib/
LPORD    = -L$(LPORDDIR) -lpord






MC30ROOT = $(current_dir)/../../thirdparty/hsl/mc30/mc30-1.0.0/src/
PARDISOROOT = /usr/local/lib/


OPENBLASROOT = /home/dav0/Apps/OpenBLAS/
BLASROOT = $(S)/BLAS/BLAS-3.7.0/

LAPACKROOT = $(S)/LAPACK/lapack-3.8.0/


LMKLIB= -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl   -m64 -I${MKLROOT}/include
LIBPARDISO = -L${PARDISOROOT} -lpardiso500-GNU481-X86-64 -lgfortran -fopenmp -lpthread -lm
# example LIBOPENBLAS= -L${OPENBLASROOT} -lopenblas_haswellp-r0.2.20.dev
LIBOPENBLAS= -L${OPENBLASROOT} -lopenblas_haswellp-r0.2.20.dev
LIBBLAS = $(BLASROOT)blas_LINUX.a
LIBLAPACK = $(LAPACKROOT)liblapack.a

LIBMC30 = ${MC30ROOT}libmc30.a -lgfortran -lm -lquadmath
K_AUG = main_mumps.c
#K_AUG = iteration4_main1.c

k_aug_pard: mc30_dri \
      get_jac_asl get_hess_asl find_ineq get_grad_f assemble_rhsds_rh assemble_rhsds_dcdp \
      suffix_decl_hand csr_driver pard_dri sig_comp mu_ad dsyev_eig dpotri_driver inertia_strategy
	$(CC) $(CCFLAGS) iteration4_main1.c -o k_aug_pardiso  $A $L \
	mc30_driver.o  \
	get_jac_asl_aug.o get_hess_asl_aug.o find_inequalities.o get_grad_f.o \
	assemble_rhs_rh.o assemble_rhs_dcdp.o  \
	suffix_decl_hand.o csr_driver.o pardiso_driver.o \
	sigma_compute.o mu_adjust_primal.o dsyev_driver.o dpotri_driver.o \
	inertia_strategy.o \
	-I$(S)/HSL/ -I$(S)/Pardiso/ -I$(S)/matrix/ -I$(S)/ASL/solvers/\
	 $(LIBMC30) $(LIBPARDISO) $(LIBLAPACK) $(LIBBLAS)

k_aug_ma57: mc30_dri \
      get_jac_asl get_hess_asl find_ineq get_grad_f assemble_rhsds_rh \
      suffix_decl_hand csr_driver sig_comp mu_ad dsyev_eig  ma57_driver
	$(CC) $(CCFLAGS) iteration4_main_ma_57.c -o k_aug_m  $A $L \
	mc30_driver.o  ma57_driver.o \
	get_jac_asl_aug.o get_hess_asl_aug.o find_inequalities.o get_grad_f.o \
	assemble_rhs_rh.o  \
	suffix_decl_hand.o csr_driver.o  \
	sigma_compute.o mu_adjust_primal.o dsyev_driver.o \
	-I$(S)/HSL/ -I$(S)/matrix/ -I$(S)/ASL/solvers/\
	 $(LIBMC30) -L/usr/local/lib/ -lcoinhsl -lmetis  $(LIBLAPACK) $(LIBBLAS)

k_aug_mumps: mc30_dri \
      get_jac_asl get_hess_asl find_ineq get_grad_f assemble_rhsds_rh assemble_rhsds_dcdp \
      suffix_decl_hand csr_driver sig_comp mu_ad dsyev_eig dpotri_driver mumps_driver inertia_strategy
	$(CC) $(CCFLAGS) $(K_AUG) -o $(ONAME)  $A $L \
	mc30_driver.o  \
	get_jac_asl_aug.o get_hess_asl_aug.o find_inequalities.o get_grad_f.o \
	assemble_rhs_rh.o assemble_rhs_dcdp.o  \
	suffix_decl_hand.o csr_driver.o mumps_driver.o \
	sigma_compute.o mu_adjust_primal.o dsyev_driver.o dpotri_driver.o \
      inertia_strategy.o \
	-I$(S)/HSL/ -I$(S)/matrix/ -I$(S)/ASL/solvers/\
	 $(LIBMC30) $(ASLINC) $(MUMPSINC) $(ASLLIB)  \
                	$(OPENBLASLIB) $(MUMPSLIB) $(LSCOTCH) $(LMETIS) $(LPORD) -lgfortran -lm -lpthread $(LIBSEQ)

k_aug_ldl: mc30_dri \
      get_jac_asl get_hess_asl find_ineq get_grad_f assemble_rhsds_rh assemble_rhsds_dcdp \
      suffix_decl_hand csr_driver pard_dri_new sig_comp mu_ad dsyev_eig dpotri_driver inertia_strategy
	$(CC) $(CCFLAGS) tstldl.c -o k_aug_ldl  $A $L \
	mc30_driver.o  \
	get_jac_asl_aug.o get_hess_asl_aug.o find_inequalities.o get_grad_f.o \
	assemble_rhs_rh.o assemble_rhs_dcdp.o  \
	suffix_decl_hand.o csr_driver.o new_pardiso_driver.o \
	sigma_compute.o mu_adjust_primal.o dsyev_driver.o dpotri_driver.o \
	inertia_strategy.o \
	-I$(S)/HSL/ -I$(S)/Pardiso/ -I$(S)/matrix/ -I$(S)/ASL/solvers/\
	 $(LIBMC30) $(LIBLAPACK) $(LIBBLAS) -ldl -lgomp



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
	$(CC) $(CCFLAGS) -c find_inequalities.c -I$S/ASL/solvers/

get_grad_f:
	$(CC) $(CCFLAGS) -c get_grad_f.c -I$S/ASL/solvers/

assemble_corrector_rhs:
	$(CC) $(CCFLAGS) -c assemble_corrector_rhs.c -I$S/ASL/solvers/

sens_update_driver:
	$(CC) $(CCFLAGS) -c $(S)/matrix/sens_update_driver.c -I$S/ASL/solvers/

assemble_rhsds_rh:
	$(CC) $(CCFLAGS) -c assemble_rhs_rh.c -I$(S)/ASL/solvers/

assemble_rhsds_dcdp:
	$(CC) $(CCFLAGS) -c assemble_rhs_dcdp.c -I$(S)/ASL/solvers/

suffix_decl_hand:
	$(CC) $(CCFLAGS) -c suffix_decl_hand.c -I$(S)/ASL/solvers/

csr_driver:
	$(CC) $(CCFLAGS) -c csr_driver.c -I$(S)/ASL/solvers/

sig_comp:
	$(CC) $(CCFLAGS) -c sigma_compute.c -I$(S)/ASL/solvers/

mu_ad:
	$(CC) $(CCFLAGS) -c mu_adjust_primal.c -I$(S)/ASL/solvers/

dsyev_eig:
	$(CC) $(CCFLAGS) -c $(S)/matrix/dsyev_driver.c -I$(S)/ASL/solvers/

ma57_driver:
	$(CC) $(CCFLAGS) -c $(S)/HSL/ma57_driver.c 
	#-L/usr/local/lib/ -lcoinhsl -lmetis $(LIBBLAS) $(LIBLAPACK) -lgfortran -lm -ldl

dpotri_driver:
	$(CC) $(CCFLAGS) -c $(S)/matrix/dpotri_driver.c 

mumps_driver:
	$(CC) $(CCFLAGS) -c ../interfaces/mumps/mumps_driver.c  $(ASLINC) $(INCSEQ)

inertia_strategy:
	$(CC) $(CCFLAGS) -c inertia_strategy.c
    
pard_dri_new:
	$(CC) $(CCFLAGS) -c $(S)/Pardiso/new_pardiso_driver.c
