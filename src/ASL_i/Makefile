CC = gcc
S = ..
L = -lm -ldl
A = $S/ASL/solvers/amplsolver.a

all1 = testin

all: $(all1)

testin: whatevs.c
	$(CC) -g  -Wall whatevs.c -o tst $A $L -I$S/ASL/solvers/
clean:
	rm -rf $(all) 

