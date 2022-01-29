#include "lu.h"
#include "mpi.h"
#include "omp.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

int main(int argc, char **argv) {
	int i, n, *perm;
	int root = 0, rank;
	double *mat, startTime, endTime;

	if (argc < 2) {
		fprintf(stderr, "Il faut donner la taille de la matrice\n");
		return 1;
	}
	n = strtol(argv[1], NULL, 10);

	srandom(time(NULL));
	mat = malloc(sizeof(double*) * n * n);
	perm = malloc(sizeof(int) * n);
	for (i = 0; i < n * n; i++) {
		mat[i] = (random() % (MAX - MIN)) + MIN;
	}

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == root) {
		startTime = MPI_Wtime();
	}
	luDecomposition(mat, n, perm);
	if (rank == root) {
		endTime = MPI_Wtime();
		printf("%lf\n", endTime - startTime);
	}
	MPI_Finalize();

	free(mat);
	free(perm);
}

void showMatrix(double *mat, int m, int n) {
	int i, j;
	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++) {
			printf("%lf ", mat[i*n + j]);
		}
		printf(";\n");
	}
}

void permute(double *A, double *tmp, int n, int i, int j) {
	if (i == j) return;
	memcpy(tmp, A + i*n, sizeof(double)*n);
	memmove(A + i*n, A + j*n, sizeof(double)*n);
	memcpy(A + j*n, tmp, sizeof(double)*n);
}

void luDecomposition(double *A, int n, int *P) {
	int i, j, k, maxIndex, startCol;
	int root = 0, size, rank, tag = 0, master;
	double maxValue, *tmp;
	MPI_Datatype column_t;
	MPI_Type_vector(n, 1, n, MPI_DOUBLE, &column_t);
	MPI_Type_commit(&column_t);

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	/* cyclic distribution of columns */
	if (rank == root) {
		j = 0;
		for (i = 0; i < n; i++) {
			if (j != rank) {
				MPI_Send(A + i, 1, column_t, j, tag, MPI_COMM_WORLD);
			}
			if (++j == size) j = 0;
		}
	} else {
		for (i = rank; i < n; i += size) {
			MPI_Recv(A + i, 1, column_t, root, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	}

	for (i = 0; i < n; i++)
		P[i] = i;

	tmp = malloc(sizeof(double) * n);
	for (i = 0; i < n; i++) {
		master = i % size;
		if (master == rank) {
			maxValue = absValue(A[i*n + i]);
			maxIndex = i;
			for (j = i; j < n; j++) {
				if (absValue(A[j*n + i]) > maxValue) {
					maxValue = absValue(A[j*n + i]);
					maxIndex = j;
				}
			}
		}
		MPI_Bcast(&maxIndex, 1, MPI_INT, master, MPI_COMM_WORLD);

		permute(A, tmp, n, i, maxIndex);
		k = P[i];
		P[i] = P[maxIndex];
		P[maxIndex] = k;

		if (master == rank) {
			//#pragma omp parallel for shared(A,n,i) private(j)
			for (j = i + 1; j < n; j++)
				A[j*n + i] /= A[i*n + i];
		}
		MPI_Bcast(A + i, 1, column_t, master, MPI_COMM_WORLD);

		startCol = i+1+(rank-master-1+size)%size;
		#pragma omp parallel for shared(A,n,i,startCol,size) private(j,k)
		for (j = i + 1; j < n; j++) {
			for (k = startCol; k < n; k += size) {
				A[j*n + k] -= A[j*n + i] * A[i*n + k];
			}
		}
	}

	/* cyclic gathering of columns */
	if (rank == root) {
		j = 0;
		for (i = 0; i < n; i++) {
			if (j != rank) {
				MPI_Recv(A + i, 1, column_t, j, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
			if (++j == size) j = 0;
		}
	} else {
		for (i = rank; i < n; i += size) {
			MPI_Send(A + i, 1, column_t, root, tag, MPI_COMM_WORLD);
		}
	}

	MPI_Type_free(&column_t);
	free(tmp);
}
