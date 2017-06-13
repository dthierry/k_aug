CC = gcc
CCFLAGS = -g -Wall -ansi -pedantic-errors -Wall -W -Wconversion -Wshadow -Wcast-qual -Wwrite-strings
FC = gfortran
S = ../../..
L = -lm -ldl
A = $S/ASL/solvers/amplsolver.a
PARENT = ../..
ONAME = dot_driver


all1 = _dot_drivers

OPENBLASROOT = /home/dav0/Apps/OpenBLAS/


LMKLIB= -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl   -m64 -I${MKLROOT}/include

# example LIBOPENBLAS= -L${OPENBLASROOT} -lopenblas_haswellp-r0.2.20.dev
LIBOPENBLAS= -L${OPENBLASROOT} -lopenblas_haswellp-r0.2.20.dev

LIBMC30 = ${MC30ROOT}libmc30.a -lgfortran -lm -lquadmath
SRC = dot_driver.c

_dot_drivers:
	$(CC) $(CCFLAGS) $(SRC) -o $(ONAME) -I$(S)/ASL/solvers/ -I$(S)/matrix/ $A $L sens_update_driver_dot.o $(LIBOPENBLAS)

sens_update_driver:
	$(CC) $(CCFLAGS) -c $(S)/matrix/sens_update_driver_dot.c -I$S/ASL/solvers/


