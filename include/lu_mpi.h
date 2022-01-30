/*! \file */
#ifndef _LU_MPI_H_
#define _LU_MPI_H_
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#define EPSILON 0.0001

/*!
 * Read square matrix from a file.
 * Reads the number of rows of the matrix on the first line.
 * Then, reads the matrix assuming that the elements are separated
 * by semicolons ';' and rows by newline character '\n'.
 *
 * \param f [in] a file opened in write mode
 * \param n [out] number of lines
 * \return matrix from the file, should be deallocated manually
 */
double *readMatrix(FILE *f, int *n);

/*!
 * Write square matrix to a file.
 * Writes the number of rows of the matrix on the first line.
 * Then, writes the matrix separating the elements by semicolons ';'
 * and rows by newline character '\n'.
 *
 * \param f [in] a file opened in write mode
 * \param mat [in] the matrix to be written
 * \param n [in] number of rows
 */
void writeMatrix(FILE *f, double *mat, int n);

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
