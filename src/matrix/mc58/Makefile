CC = gcc
CCFLAGS = -g -Wall -ansi -pedantic-errors
SRC = mc58driver.c
L = -lm -ldl
MC58 = /home/dav0/Apps/mc58-1.0.0/src/libmc58.a
S = ../..
LBLAS =  $(S)/BLAS/BLAS-3.7.0/blas_LINUX.a

all: mc58_driver

mc58_driver:
	$(CC) $(CCFLAGS) $(SRC) -o mc58_driver_ $(MC58) $(LBLAS) -lgfortran



