CC = gcc
FC = gfortran
S = ..
L = -lm -ldl
A = $S/ASL/solvers/amplsolver.a
ONAME = k_matrix


all1 = k_matrix_openblas
MATRIXDIR = ../matrix
LMATRIX = kassemble.f

# example MC30ROOT = /home/dav0/Apps/mc30-1.0.0/build/src/
MC30ROOT = 

# example PARDISOROOT = /usr/local/lib/
PARDISOROOT = 

# example OPENBLASROOT = /home/dav0/Apps/OpenBLAS/
OPENBLASROOT = 

LMKLIB= -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl   -m64 -I${MKLROOT}/include
LIBPARDISO = -L${PARDISOROOT} -lpardiso500-GNU481-X86-64 -lgfortran -fopenmp -lpthread -lm
# example LIBOPENBLAS= -L${OPENBLASROOT} -lopenblas_haswellp-r0.2.20.dev
LIBOPENBLAS= -L${OPENBLASROOT} -lopenblas_haswellp-r0.2.20.dev

LIBMC30 = ${MC30ROOT}libmc30.a -lgfortran -lm -lquadmath
KMATRIX = kmatrix_corrector.c



all: $(all1)
k_matrix_mlk: 
	$(CC) -g  -Wall $(KMATRIX) -o test1  $A $L \
	-I$S/ASL/solvers/ \
	sort_by_column.o w_append_nz.o k_assemble_cc.o kmalloc.o mc30_driver.o  \
	get_jac_asl.o get_hess_asl.o find_inequalities.o get_grad_f.o \
	assemble_corrector_rhs.o assemble_rhsds.o \
	$(MC30ROOT)libmc30.a \
      -lgfortran -lm -lquadmath \
      -I$(S)/HSL/ -I$(S)/Pardiso/ $(LIBPARDISO) $(LMKLIB) $(S)/Pardiso/pardiso_driver.o
k_matrix_openblas: sort_by_column wzndiag k_assemble kmalloc mc30_dri pard_dri \
      get_jac_asl get_hess_asl find_ineq get_grad_f assemble_corrector_rhs assemble_rhsds_v3
	$(CC) -g  -Wall $(KMATRIX) -o $(ONAME)  $A $L \
	-I$S/ASL/solvers/ \
	sort_by_column.o w_append_nz.o k_assemble_cc.o kmalloc.o mc30_driver.o  \
	get_jac_asl.o get_hess_asl.o find_inequalities.o get_grad_f.o \
	assemble_corrector_rhs.o assemble_rhsds.o pardiso_driver.o \
	-I$(S)/HSL/ -I$(S)/Pardiso/ \
	$(LIBPARDISO) $(LIBOPENBLAS) $(LIBMC30) 


clean:
	rm -rf $(all) *.o
sort_by_column:
	$(CC) -g -Wall -c sort_by_column.c -I$S/ASL/solvers/
wzndiag: 
	$(CC) -g -Wall -c w_append_nz.c -I$S/ASL/solvers/
k_assemble:
	$(CC) -g -Wall -c k_assemble_cc.c -I$S/ASL/solvers/ $L
kmalloc:
	$(CC) -g -Wall -c kmalloc.c -I$S/ASL/solvers/
      
mc30_dri:
	$(CC) -g -Wall -c ../HSL/mc30_driver.c -I$(S)/HSL/ -I$S/ASL/solvers/

pard_dri:
	$(CC) -g -Wall -c $(S)/Pardiso/pardiso_driver.c -I$(PWD)/../ASL/solvers/ 


get_jac_asl:
	$(CC) -g -Wall -c get_jac_asl.c -I$S/ASL/solvers/

get_hess_asl:
	$(CC) -g -Wall -c get_hess_asl.c -I$S/ASL/solvers/

find_ineq:
	$(CC) -g -Wall -c find_inequalities.c -I$S/ASL/solvers/


get_grad_f:
	$(CC) -g -Wall -c get_grad_f.c -I$S/ASL/solvers/

assemble_corrector_rhs:
	$(CC) -g -Wall -c assemble_corrector_rhs.c -I$S/ASL/solvers/
assemble_rhsds:
	$(CC) -g -Wall -c assemble_rhsds.c -o assemble_rhsds.o -I$(S)/ASL/solvers/
assemble_rhsds_v2:
	$(CC) -g -Wall -c assemble_rhsdsv2.c -o assemble_rhsds.o -I$(S)/ASL/solvers/
assemble_rhsds_v3:
	$(CC) -g -Wall -c assemble_rhsdsv3.c -o assemble_rhsds.o -I$(S)/ASL/solvers/
