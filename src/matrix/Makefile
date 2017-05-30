CC = gcc
CCFLAGS = -g -Wall -ansi -pedantic-errors
SRC = vec_mat_prd.ca
S = ..
L = -lm -ldl
A = $S/ASL/solvers/amplsolver.a

OPENBLASROOT = /home/dav0/Apps/OpenBLAS/
LIBOPENBLAS= -L${OPENBLASROOT} -lopenblas_haswellp-r0.2.20.dev

vm:
	$(CC) $(CCFLAGS) vec_mat_prd.c $(LIBOPENBLAS) -lgfortran -ldl -lm

driver:
	$(CC) $(CCFLAGS) -c sens_update_driver.c -I$S/ASL/solvers/
