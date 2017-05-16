CC = gcc
CCFLAGS = -g -Wall -ansi -pedantic-errors
SRC = vec_mat_prd.ca

OPENBLASROOT = /home/dav0/Apps/OpenBLAS/
LIBOPENBLAS= -L${OPENBLASROOT} -lopenblas_haswellp-r0.2.20.dev

vm:
	$(CC) $(CCFLAGS) vec_mat_prd.c $(LIBOPENBLAS) -lgfortran -ldl -lm
