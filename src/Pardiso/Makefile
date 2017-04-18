FC=gfortran
FCFLAGS=-O2 -g
CC=gcc
CCFLAGS=-O2 -g -c

SRCFILE=pardiso_driver.c
# threaded layer 32 bit # works

LMKLIB= -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -lgfortran -lgomp -lpthread -lm -ldl -m64 -I${MKLROOT}/include

LIBPARDISO = -L/usr/local/lib/ -lpardiso500-GNU481-X86-64

all: pard_test 
pard_test:
	$(CC) $(CCFLAGS) $(SRCFILE) $(LIBPARDISO) $(LMKLIB) -I$(PWD)/../ASL/solvers/

