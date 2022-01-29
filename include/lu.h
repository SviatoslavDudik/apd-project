/*! \file */
#ifndef _LU_H_
#define _LU_H_

/*! Minimum value to generate in the matrix */
#define MIN -10
/*! Maximum value to generate in the matrix */
#define MAX 10
/*! Macro to compute the absolute value of a number */
#define absValue(x) ((x) >= 0 ? (x) : -(x))

/*!
 * Print matrix.
 *
 * Print matrix separating the elements by spaces and lines by semicolons ';'.
 *
 * \param mat [in] the matrix to be outputed
 * \param m [in] number of lines
 * \param n [in] number of columns
 */
void showMatrix(double *mat, int m, int n);

/*!
 * Permute two lines of a matrix.
 *
 * Exchanges line i with line j in a square matrix of size n.
 *
 * \param A [in,out] matrix
 * \param tmp [in] temporary array of size n
 * \param n matrix size
 * \param i [in] first line number
 * \param j [in] second line number
 */
void permute(double *A, double *tmp, int n, int i, int j);

/*!
 * Performs LU factorization with partial pivoting.
 *
 * Decomposes matrix A in 2 matrices L and U so that
 * PA = LU, where:
 * L - lower diagonal matrix with ones on the main diagonal;
 * U - upper diagonal matrix;
 * P - permutation matrix.
 *
 * A must be an array of size n*n.
 * P must be an array of size n.
 *
 * At the end of the function, A will be equal to U + (L - I) and P will contain
 * values from 0 to (n-1) giving the permutation to be performed on the initial
 * matrix.
 *
 * \param A [in,out] input matrix and also the output matrix containing U + (L - I)
 * \param n [in] number of lines
 * \param P [out] permutation vector
 */
void luDecomposition(double *A, int n, int *P);

#endif
