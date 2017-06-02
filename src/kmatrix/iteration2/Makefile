CC = gcc
CCFLAGS = -g -Wall -ansi -pedantic-errors -Wall -W -Wconversion -Wshadow -Wcast-qual -Wwrite-strings
FC = gfortran
S = ../..
L = -lm -ldl
A = $S/ASL/solvers/amplsolver.a
PARENT = ..
ONAME = k_aug


all1 = k_aug
MATRIXDIR = ../matrix
LMATRIX = kassemble.f

MC30ROOT = /home/dav0/Apps/mc30-1.0.0/build/src/


PARDISOROOT = /usr/local/lib/


OPENBLASROOT = /home/dav0/Apps/OpenBLAS/


LMKLIB= -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl   -m64 -I${MKLROOT}/include
LIBPARDISO = -L${PARDISOROOT} -lpardiso500-GNU481-X86-64 -lgfortran -fopenmp -lpthread -lm
# example LIBOPENBLAS= -L${OPENBLASROOT} -lopenblas_haswellp-r0.2.20.dev
LIBOPENBLAS= -L${OPENBLASROOT} -lopenblas_haswellp-r0.2.20.dev

LIBMC30 = ${MC30ROOT}libmc30.a -lgfortran -lm -lquadmath
KMATRIX = iteration1_main.c


k_aug: mc30_dri \
      get_jac_asl get_hess_asl find_ineq get_grad_f assemble_rhsds_rh \
      suffix_decl_hand csr_driver pard_dri
	$(CC) $(CCFLAGS) iteration2_main1.c -o $(ONAME)  $A $L \
	mc30_driver.o  \
	get_jac_asl_aug.o get_hess_asl_aug.o find_inequalities.o get_grad_f.o \
	assemble_rhsd_red_hess.o  \
	suffix_decl_hand.o csr_driver.o pardiso_driver.o \
	-I$(S)/HSL/ -I$(S)/Pardiso/ -I$(S)/matrix/ -I$(PARENT) -I$(S)/ASL/solvers/\
	 $(LIBMC30) $(LIBPARDISO) $(LMKLIB)

k_aug_oblas: mc30_dri \
      get_jac_asl get_hess_asl find_ineq get_grad_f assemble_rhsds_rh \
      suffix_decl_hand csr_driver pard_dri
	$(CC) $(CCFLAGS) iteration2_main1.c -o $(ONAME)  $A $L \
	mc30_driver.o  \
	get_jac_asl_aug.o get_hess_asl_aug.o find_inequalities.o get_grad_f.o \
	assemble_rhsd_red_hess.o  \
	suffix_decl_hand.o csr_driver.o pardiso_driver.o \
	-I$(S)/HSL/ -I$(S)/Pardiso/ -I$(S)/matrix/ -I$(PARENT) -I$(S)/ASL/solvers/\
	 $(LIBMC30) $(LIBPARDISO) $(LIBOPENBLAS)



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
	$(CC) $(CCFLAGS) -c $(PARENT)/assemble_corrector_rhs.c -I$S/ASL/solvers/

assemble_rhsds_v3:
	$(CC) $(CCFLAGS) -c $(PARENT)/assemble_rhsdsv3.c -o assemble_rhsds.o -I$(S)/ASL/solvers/

sens_update_driver:
	$(CC) $(CCFLAGS) -c $(S)/matrix/sens_update_driver.c -I$S/ASL/solvers/

assemble_rhsds_rh:
	$(CC) $(CCFLAGS) -c $(PARENT)/assemble_rhsd_red_hess.c -I$(S)/ASL/solvers/

suffix_decl_hand:
	$(CC) $(CCFLAGS) -c suffix_decl_hand.c -I$(S)/ASL/solvers/

csr_driver:
	$(CC) $(CCFLAGS) -c csr_driver.c -I$(S)/ASL/solvers/