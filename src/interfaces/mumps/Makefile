FC=gfortran
FCFLAGS= -O2 -g
CC=gcc
CCFLAGS= -O2 -g
DEPS =
OPENBLASLIB = -L../../../thirdparty/openblas/OpenBLAS/ -lopenblas
MUMPSLIB = -L/home/dav0/k2/KKT_matrix/thirdparty/mumps512/MUMPS_5.1.2/lib/ -ldmumps -lmumps_common
MUMPSINC = -I/home/dav0/k2/KKT_matrix/thirdparty/mumps512/MUMPS_5.1.2/include/
ASLLIB = ../../../src/ASL/solvers/amplsolver.a
ASLINC = -I/home/dav0/k2/KKT_matrix/src/ASL/solvers/
LIBSEQ = -L/home/dav0/k2/KKT_matrix/thirdparty/mumps512/MUMPS_5.1.2/libseq/ -lmpiseq
INCSEQ = -I/home/dav0/k2/KKT_matrix/thirdparty/mumps512/MUMPS_5.1.2/libseq/

LSCOTCHDIR = /home/dav0/k2/KKT_matrix/thirdparty/scotch/whatnot/scotch_6.0.4/lib/
LSCOTCH   = -L$(LSCOTCHDIR) -lesmumps -lscotch -lscotcherr
LMETISDIR = /home/dav0/k2/KKT_matrix/thirdparty/metis/metis-5.1.0/build/Linux-x86_64/libmetis/
LMETIS    = -L$(LMETISDIR) -lmetis

LPORDDIR = -L/home/dav0/k2/KKT_matrix/thirdparty/mumps512/MUMPS_5.1.2/PORD/lib/
LPORD    = -L$(LPORDDIR) -lpord



all: mumps_test

mumps_driver:
	$(CC) $(CCFLAGS) -c mumps_driver.c  $(ASLINC) $(INCSEQ)

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $<

mumps_test: mumps_driver
	$(CC) $(CCFLAGS) -o whatevs test_mumps.c mumps_driver.o $(ASLINC) $(MUMPSINC) $(ASLLIB)  \
	$(OPENBLASLIB) $(MUMPSLIB) $(LSCOTCH) $(LMETIS) $(LPORD) -lgfortran -lm -lpthread $(LIBSEQ)

.PHONY: clean

clean:
	rm -f *.o
