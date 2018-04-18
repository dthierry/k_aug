FC=gfortran
FCFLAGS= -O2 -g
CC=gcc
CCFLAGS= -O2 -g
DEPS =
OPENBLASLIB = /usr/local/lib/openblas -lopenblas
MUMPSLIB = -L/home/dav0/k2/KKT_matrix/thirdparty/mumps512/MUMPS_5.1.2/lib -ldmumps
MUMPSINC = -I/home/dav0/k2/KKT_matrix/thirdparty/mumps512/MUMPS_5.1.2/include/
ASLINC = -I/home/dav0/k2/KKT_matrix/src/ASL/solvers/

all: mumps_test

mumps_driver:
	$(CC) $(CCFLAGS) -c mumps_driver.c  $(ASLINC)

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $<

mumps_test: mumps_driver
	$(CC) $(CCFLAGS) -o whatevs test_mumps.c mumps_driver.o  $(ASLINC) $(MUMPSINC) ../../ASL/solvers/amplsolver.a

.PHONY: clean

clean:
	rm -f *.o
