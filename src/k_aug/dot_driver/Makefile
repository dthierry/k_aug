CC = gcc
CCFLAGS = -g -Wall -ansi -pedantic-errors -Wall -W -Wconversion -Wshadow -Wcast-qual -Wwrite-strings -ggdb
FC = gfortran
S = ../..
L = -lm -ldl -lgfortran
A = -L$(current_dir)/../../../thirdparty/asl/solvers/ -lamplsolver

ONAME = dot_driver
all1 = _dot_drivers

mkfile_path := $(abspath $(lastword $(MAKEFILE_LIST)))
current_dir := $(patsubst %/,%,$(dir $(mkfile_path)))

OPENBLASLIB = -L$(current_dir)/../../../thirdparty/openblas/OpenBLAS/ -lopenblas

LIBMC30 = ${MC30ROOT}libmc30.a -lgfortran -lm -lquadmath
SRC = dot_driver.c

_dot_drivers:
	$(CC) $(CCFLAGS) $(SRC) -o $(ONAME) -I$(S)/ASL/solvers/ -I$(S)/matrix/ $A $(OPENBLASLIB) $L


