FC=gfortran
FCFLAGS= -O2 -g
CC=gcc
CCFLAGS= -O2 -g
DEPS = mumps_driver.h
OPENBLASLIB = /usr/local/lib/openblas -lopenblas
MUMPSLIB = /home/dav0/Apps/MUMPS_5.1.2/lib -ldmumps
MUMPSINC = -I/home/dav0/Apps/MUMPS_5.1.2/include
ASLINC = -I../../ASL/solvers

all: mumps_test

#mumps_driver:
#	$(CC) $(CCFLAGS) -c mumps_driver.c $(ASLINC) $(MUMPSINC)

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CCFLAGS) $(ASLINC)

mumps_test: mumps_driver.o
	$(CC) $(CCFLAGS) -o whatevs test_mumps.c mumps_driver.o $(ASLINC) $(MUMPSINC) -I.

.PHONY: clean

clean:
	rm -f *.o
