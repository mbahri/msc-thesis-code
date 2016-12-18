# makefile
# Makefile for compiling shrinkage_p.cpp.
#
# Georgios Papamakarios
# Imperial College London
# Aug 2014

CC = gcc
MATINC = /usr/lib/matlab/R2013a/extern/include
CFLAGS = -fPIC -O3

shrinkage_p.o: shrinkage_p.cpp
	$(CC) -c shrinkage_p.cpp -o shrinkage_p.o $(CFLAGS) -I$(MATINC)

clean:
	rm -f shrinkage_p.o
