CC = gcc
S = ..
L = -lm -ldl
A = $S/ASL/solvers/amplsolver.a

all1 = testin

all: $(all1)

testin: kmatrix.c
	$(CC) -g  -Wall kmatrix.c -o tst $A $L -I$S/ASL/solvers/
clean:
	rm -rf $(all) 

