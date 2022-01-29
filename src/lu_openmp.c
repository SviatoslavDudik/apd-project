#include "lu_openmp.h"
#include "omp.h"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) {
	int i, j, n, *perm;
	double **mat, startTime, endTime;

	if (argc < 2) {
		fprintf(stderr, "Il faut donner la taille de la matrice\n");
		return 1;
	}
	n = strtol(argv[1], NULL, 10);

	mat = malloc(sizeof(double*) * n);
	perm = malloc(sizeof(int) * n);
	for (i = 0; i < n; i++) {
		mat[i] = malloc(sizeof(double) * n);
		for (j = 0; j < n; j++) {
			mat[i][j] = (random() % (1000*(MAX - MIN)))/1000 + MIN;
		}
	}

	startTime = omp_get_wtime();
	luDecomposition(mat, n, perm);
	endTime = omp_get_wtime();

	printf("%lf\n", endTime - startTime);
	for (i = 0; i < n; i++)
		free(mat[i]);
	free(mat);
}

void showMatrix(double **mat, int m, int n) {
	int i, j;
	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++) {
			printf("%lf ", mat[i][j]);
		}
		printf("\n");
	}
}

void permute(double **A, int i, int j) {
	double *temp;
	temp = A[i];
	A[i] = A[j];
	A[j] = temp;
}

void luDecomposition(double **A, int n, int *P) {
	int i, j, k, maxIndex;
	double maxValue;
	for (i = 0; i < n; i++)
		P[i] = i;
	for (i = 0; i < n; i++) {
		maxValue = absValue(A[i][i]);
		maxIndex = i;
		for (j = i; j < n; j++) {
			if (absValue(A[j][i]) > maxValue) {
				maxValue = absValue(A[j][i]);
				maxIndex = j;
			}
		}
		permute(A, i, maxIndex);
		k = P[i];
		P[i] = P[maxIndex];
		P[maxIndex] = k;

		#pragma omp parallel for shared(A,n,i) private(j,k)
		for (j = i + 1; j < n; j++) {
			A[j][i] /= A[i][i];
			for (k = i + 1; k < n; k++) {
				A[j][k] -= A[j][i] * A[i][k];
			}
		}
	}
}

