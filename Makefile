MPICC=mpicc
CC=gcc
CFLAGS=-Wall -O3
OMPFLAGS=-fopenmp

.PHONY=all

executables = lu lu_mpi lu_openmp monte monte_mpi monte_openmp
all: $(executables)
clean:
	rm $(executables)

lu: lu.c lu.h
	$(MPICC) $(CFLAGS) $(OMPFLAGS) $< -o $@

lu_mpi: lu_mpi.c lu_mpi.h
	$(MPICC) $(CFLAGS) $< -o $@

lu_openmp: lu_openmp.c lu_openmp.h
	$(CC) $(CFLAGS) $(OMPFLAGS) $< -o $@

monte: monte.c
	$(MPICC) $(CFLAGS) $(OMPFLAGS) $< -o $@

monte_mpi: monte_mpi.c
	$(MPICC) $(CFLAGS) $< -o $@

monte_openmp: monte_openmp.c
	$(CC) $(CFLAGS) $(OMPFLAGS) $< -o $@
