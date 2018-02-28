CC = gcc
CCFLAGS = -g -Wall -ansi -pedantic-errors -Wall -W -Wconversion -Wshadow -Wcast-qual -Wwrite-strings
FC = gfortran
S = ../..
L = -lm -ldl -lgfortran
A = $S/ASL/solvers/amplsolver.a

ONAME = dot_driver
all1 = _dot_drivers

OPENBLASROOT = /home/dav0/Apps/OpenBLAS/

BLASROOT = $(S)/BLAS/BLAS-3.7.0/
LIBBLAS = $(BLASROOT)blas_LINUX.a

LIBMC30 = ${MC30ROOT}libmc30.a -lgfortran -lm -lquadmath
SRC = dot_driver.c

_dot_drivers: sens_update_driver
	$(CC) $(CCFLAGS) $(SRC) -o $(ONAME) -I$(S)/ASL/solvers/ -I$(S)/matrix/ $A sens_update_driver_dot.o $(LIBBLAS) $L

sens_update_driver:
	$(CC) $(CCFLAGS) -c $(S)/matrix/sens_update_driver_dot.c -I$(S)/ASL/solvers/ 


