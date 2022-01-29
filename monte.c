#include "mpi.h"
#include "omp.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define RADIUS 1000

int main(int argc, char **argv) {
	int root = 0, rank, size;
	long n, i, count = 0, countTotal;
	double x, y, pi, startTime, endTime;
	n = strtol(argv[1], NULL, 10);

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	srandom(time(NULL));
	MPI_Barrier(MPI_COMM_WORLD);
	if (rank == root) {
		startTime = MPI_Wtime();
	}
	#pragma omp parallel for private(i,x,y) shared(n,size) reduction(+:count)
	for (i = 0; i < 1000 * n / size; i++) {
		x = (double)(random() % (RADIUS*1000))/1000;
		y = (double)(random() % (RADIUS*1000))/1000;
		if (x*x + y*y < RADIUS*RADIUS) {
			count++;
		}
	}
	MPI_Reduce(&count, &countTotal, 1, MPI_LONG, MPI_SUM, root, MPI_COMM_WORLD);

	if (rank == root) {
		pi = (double) 4 * countTotal / (1000 * n);
		endTime = MPI_Wtime();
		printf("%lf\n", endTime - startTime);
	}
	MPI_Finalize();
	return EXIT_SUCCESS;
}
