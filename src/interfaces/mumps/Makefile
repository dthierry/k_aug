FC=gfortran
FCFLAGS= -O2 -g
CC=gcc
CCFLAGS= -O2 -g

OPENBLASLIB = /usr/local/lib/openblas -lopenblas
MUMPSLIB = /home/dav0/Apps/MUMPS_5.1.2/lib -ldmumps
MUMPSINC = -I/home/dav0/Apps/MUMPS_5.1.2/include
ASLINC = -I/home/dav0/Apps/k_aug/src/ASL/solvers

all: mumps_driver mumps_test

mumps_driver:
	$(CC) $(CCFLAGS) -c mumps_driver.c $(ASLINC) $(MUMPSINC)

mumps_test:
	$(CC) $(CCFLAGS) -o whatevs test_mumps.c mumps_driver.o $(ASLINC) $(MUMPSINC)