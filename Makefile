MPICC=mpicc
CC=gcc
CFLAGS=-Wall -O3 -lm
OMPFLAGS=-fopenmp

SRCDIR=src
INCLUDEDIR=include

executables = lu lu_mpi lu_openmp monte monte_mpi monte_openmp

.PHONY=all

all: $(executables)

clean:
	rm $(executables)

lu: $(SRCDIR)/lu.c $(INCLUDEDIR)/lu.h
	$(MPICC) $(CFLAGS) $(OMPFLAGS) -I $(INCLUDEDIR) $< -o $@

lu_mpi: $(SRCDIR)/lu_mpi.c $(INCLUDEDIR)/lu_mpi.h
	$(MPICC) $(CFLAGS) -I $(INCLUDEDIR) $< -o $@

lu_openmp: $(SRCDIR)/lu_openmp.c $(INCLUDEDIR)/lu_openmp.h
	$(CC) $(CFLAGS) $(OMPFLAGS) -I $(INCLUDEDIR) $< -o $@

monte: $(SRCDIR)/monte.c
	$(MPICC) $(CFLAGS) $(OMPFLAGS) $< -o $@

monte_mpi: $(SRCDIR)/monte_mpi.c
	$(MPICC) $(CFLAGS) $< -o $@

monte_openmp: $(SRCDIR)/monte_openmp.c
	$(CC) $(CFLAGS) $(OMPFLAGS) $< -o $@
