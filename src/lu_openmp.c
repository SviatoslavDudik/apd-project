#include "lu_openmp.h"

int main(int argc, char **argv) {
	FILE *f;
	int i, n, *perm;
	double **mat, startTime, endTime;

	if (argc < 2) {
		fprintf(stderr, "Il faut donner le fichier avec la matrice\n");
		return 1;
	}
	f = fopen(argv[1], "r");
	if (f == NULL) {
		fprintf(stderr, "Impossible de lire le fichier\n");
		return 1;
	}
	mat = readMatrix(f, &n);
	fclose(f);
	perm = malloc(sizeof(int) * n);

	startTime = omp_get_wtime();
	luDecomposition(mat, n, perm);
	endTime = omp_get_wtime();

	printf("%lf\n", endTime - startTime);
	for (i = 0; i < n; i++)
		free(mat[i]);
	free(mat);
}

double **readMatrix(FILE *f, int *n) {
	int i, j;
	double **mat;
	fscanf(f, "%d\n", n);
	mat = malloc(sizeof(double*) * (*n));
	for (i = 0; i < *n; i++) {
		mat[i] = malloc(sizeof(double) * (*n));
		for (j = 0; j < *n; j++) {
			fscanf(f, "%lf;", mat[i] + j);
		}
		fscanf(f, "\n");
	}
	return mat;
}

void writeMatrix(FILE *f, double **mat, int n) {
	int i, j;
	fprintf(f, "%d\n", n);
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			fprintf(f, "%lf;", mat[i][j]);
		}
		fprintf(f, "\n");
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
		maxValue = fabs(A[i][i]);
		maxIndex = i;
		for (j = i; j < n; j++) {
			if (fabs(A[j][i]) > maxValue) {
				maxValue = fabs(A[j][i]);
				maxIndex = j;
			}
		}
		if (maxValue < EPSILON)
			fprintf(stderr, "ERROR: matrix is degenerate\n");
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

