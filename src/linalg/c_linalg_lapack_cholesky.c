/*
 * c_linalg.c
 *
 *  Created on: 2014/05/13
 *      Author: utsugi
 */

#include <math.h>
#include <clinalg.h>

#include "private.h"

/* interface of lapack dpotrf_ */
int
c_linalg_lapack_dpotrf (char uplo, c_matrix *a)
{
	int		info;
	int		n;
	int		lda;

	if (c_matrix_is_empty (a)) c_error ("c_linalg_lapack_dpotrf", "matrix is empty.");
	if (!c_matrix_is_square (a)) c_error ("c_linalg_lapack_dpotrf", "matrix must be square.");
	if (uplo != 'L' && uplo != 'U') c_error ("c_linalg_lapack_dpotrf", "uplo must be 'L' or 'U'.");

	n = (int) a->size1;
	lda = (int) a->lda;
	F77CALL (dpotrf) (&uplo, &n, a->data, &lda, &info);
	return info;
}

/* interface of lapack dpotrs_ */
int
c_linalg_lapack_dpotrs (char uplo, c_matrix *l, c_matrix *b)
{
	int		info;
	int		n;
	int		nrhs;
	int		lda;
	int		ldb;

	if (c_matrix_is_empty (l)) c_error ("c_linalg_lapack_dpotrs", "first matrix is empty.");
	if (c_matrix_is_empty (b)) c_error ("c_linalg_lapack_dpotrs", "second matrix is empty.");
	if (!c_matrix_is_square (l)) c_error ("c_linalg_lapack_dpotrs", "first matrix must be square.");
	if (l->size1 != b->size1) c_error ("c_linalg_lapack_dpotrs", "matrix size does not match.");
	if (uplo != 'L' && uplo != 'U') c_error ("c_linalg_lapack_dpotrs", "uplo must be 'L' or 'U'.");

	n = (int) l->size1;
	nrhs = (int) b->size2;
	lda = (int) l->lda;
	ldb = (int) b->lda;
	F77CALL (dpotrs) (&uplo, &n, &nrhs, l->data, &lda, b->data, &ldb, &info);
	return info;
}

/* interface of lapack dpotri_ */
int
c_linalg_lapack_dpotri (char uplo, c_matrix *l)
{
	int		info;
	int		n;
	int		lda;
	if (c_matrix_is_empty (l)) c_error ("c_linalg_lapack_dpotri", "matrix is empty.");
	if (!c_matrix_is_square (l)) c_error ("c_linalg_lapack_dpotri", "matrix must be square.");
	if (uplo != 'L' && uplo != 'U') c_error ("c_linalg_lapack_dpotri", "uplo must be 'L' or 'U'.");

	n = (int) l->size1;
	lda = (int)l->lda;
	F77CALL (dpotri) (&uplo, &n, l->data, &lda, &info);
	return info;
}

