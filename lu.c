#include <stdio.h>
#include <stdlib.h>

#define abs(x) ((x) >= 0 ? (x) : -(x))

void showMatrix(double **mat, int m, int n);
void permute(double **A, int n, int i, int j);
void luDecomposition(double **A, int n, int *P);

int main(int argc, char **argv) {
	int i, j, n, *perm;
	FILE *f;
	double **mat;

	if (argc < 2) {
		fprintf(stderr, "Il faut donner le nom de fichier avec la matrice\n");
		return 1;
	}
	f = fopen(argv[1], "r");
	if (f == NULL) {
		fprintf(stderr, "Impossible d'ouvrir le fichier\n");
		return 1;
	}
	fscanf(f, "%d", &n);

	mat = malloc(sizeof(double*) * n);
	perm = malloc(sizeof(int) * n);
	for (i = 0; i < n; i++) {
		mat[i] = malloc(sizeof(double) * n);
		for (j = 0; j < n; j++) {
			fscanf(f, "%lf", mat[i] + j);
		}
	}

	showMatrix(mat, n, n);
	luDecomposition(mat, n, perm);
	showMatrix(mat, n, n);

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

void permute(double **A, int n, int i, int j) {
	int k;
	double temp;
	for (k = 0; k < n; k++) {
		temp = A[i][k];
		A[i][k] = A[j][k];
		A[j][k] = temp;
	}
}

void luDecomposition(double **A, int n, int *P) {
	int i, j, k, maxIndex;
	double maxValue;
	for (i = 0; i < n; i++) {
		maxValue = abs(A[i][i]);
		maxIndex = i;
		for (j = i; j < n; j++) {
			if (abs(A[j][i]) > maxValue) {
				maxValue = abs(A[j][i]);
				maxIndex = j;
			}
		}
		permute(A, n, i, maxIndex);
		P[i] = maxIndex;

		for (j = i + 1; j < n; j++) {
			A[j][i] /= A[i][i];
			for (k = i + 1; k < n; k++) {
				A[j][k] -= A[j][i] * A[i][k];
			}
		}
	}
}

