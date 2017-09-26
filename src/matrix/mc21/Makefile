CC = gcc
CCFLAGS = -g -Wall -ansi -pedantic-errors
SRC = mc21.c
S = ../..
MC21 = /home/dav0/Apps/mc21-1.0.0/build/src/libmc21.a

all: mc21_driver

mc21_driver:
	$(CC) $(CCFLAGS) $(SRC) -o mc21_driver_ $(MC21) -lgfortran
