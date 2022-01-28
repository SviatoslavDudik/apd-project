#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MIN -10
#define MAX 10
#define abs(x) ((x) >= 0 ? (x) : -(x))

void showMatrix(double *mat, int m, int n);
void memfill(int *arr, int n, int value);
void permute(double *A, double *tmp, int n, int i, int j);
void luDecomposition(double *A, int n, int *P);

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
		printf("Time: %lf s\n", endTime - startTime);
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
		printf("\n");
	}
}

void permute(double *A, double *tmp, int n, int i, int j) {
	if (i == j) return;
	memcpy(tmp, A + i*n, sizeof(double)*n);
	memmove(A + i*n, A + j*n, sizeof(double)*n);
	memcpy(A + j*n, tmp, sizeof(double)*n);
}

void memfill(int *arr, int n, int value) {
	int *arrdest = arr + n;
	while (arr < arrdest)
		*arr++ = value;
}

void luDecomposition(double *A, int n, int *P) {
	int i, j, k, maxIndex, row;
	int root = 0, size, rank;
	double maxValue, *tmp, *start;
	int *counts, *displs, rows, chunk;
	MPI_Datatype ROW;
	MPI_Type_contiguous(n, MPI_DOUBLE, &ROW);
	MPI_Type_commit(&ROW);

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	tmp = malloc(sizeof(double) * n);
	counts = malloc(sizeof(int) * size);
	displs = malloc(sizeof(int) * size);
	for (i = 0; i < n; i++) {
		if (rank == root) {
			maxValue = abs(A[i*n + i]);
			maxIndex = i;
			for (j = i; j < n; j++) {
				if (abs(A[j*n + i]) > maxValue) {
					maxValue = abs(A[j*n + i]);
					maxIndex = j;
				}
			}
			permute(A, tmp, n, i, maxIndex);
			P[i] = maxIndex;
		}

		MPI_Bcast(A + i*n + i, n - i, MPI_DOUBLE, root, MPI_COMM_WORLD);
		rows = n - i - 1;
		chunk = rows / (size - 1);
		memfill(counts, size, chunk);
		counts[root] = rows % (size - 1);
		displs[0] = 0;
		for (j = 1; j < size; j++) {
			displs[j] = displs[j-1] + counts[j-1];
		}
		start = A + (i+1) * n;

		if (rank == root)
			MPI_Scatterv(start, counts, displs, ROW, MPI_IN_PLACE, counts[rank], ROW, root, MPI_COMM_WORLD);
		else
			MPI_Scatterv(MPI_IN_PLACE, counts, displs, ROW, start, counts[rank], ROW, root, MPI_COMM_WORLD);
		for (j = 0; j < counts[rank]; j++) {
			row = (i + 1) + displs[rank] + j;
			A[row*n + i] /= A[i*n + i];
			for (k = i + 1; k < n; k++) {
				A[row*n + k] -= A[row*n + i] * A[i*n + k];
			}
		}
		if (rank == root)
			MPI_Gatherv(MPI_IN_PLACE, counts[rank], ROW, start, counts, displs, ROW, root, MPI_COMM_WORLD);
		else
			MPI_Gatherv(start, counts[rank], ROW, MPI_IN_PLACE, counts, displs, ROW, root, MPI_COMM_WORLD);
	}
	MPI_Type_free(&ROW);
	free(tmp);
	free(counts);
	free(displs);
}
