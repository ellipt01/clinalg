/*
 * c_matrixops.c
 *
 *  Created on: 2014/04/03
 *      Author: utsugi
 */

#include <c_matrix.h>

extern void	c_error (const char * function_name, const char *error_msg);

/* blas */
extern void	dcopy_ (long *n, double *x, long *incx, double *y, long *incy);
extern void	dgemv_ (char *trans, long *n, long *m, double *alpha, double *a, long *lda, double *x, long *incx, double *beta, double *y, long *incy);
extern void	dgemm_ (char *transA, char *transB, long *m, long *n, long *k, double *alpha, double *a, long *lda, double *b, long *ldb, double *beta, double *c, long *ldc);

c_matrix *
c_matrix_transpose (c_matrix *a)
{
	int			j;
	c_vector	*col;
	c_matrix	*at = c_matrix_alloc (a->size2, a->size1);

	long		n;
	long		incx;
	long		incy = (long) at->lda;
	for (j = 0; j < a->size2; j++) {
		col = c_matrix_column (a, j);
		n = (long) a->size1;
		incx = (long) col->stride;
		dcopy_ (&n, col->data, &incx, at->data + j, &incy);
	}
	return at;
}

/* y = alpha * a * x + beta */
c_vector *
c_matrix_dot_vector (double alpha, const c_matrix *a, const c_vector *x, double beta)
{
	char		trans = 'N';
	long		n;
	long		m;
	long		lda;
	long		incx;
	long		incy;
	c_vector	*y;
	if (c_matrix_is_empty (a)) c_error ("c_matrix_dot_vector", "matrix is empty.");
	if (c_vector_is_empty (x)) c_error ("c_matrix_dot_vector", "vector is empty.");
	if (a->size2 != x->size) c_error ("c_matrix_dot_vector", "vector and matrix size does not match.");
	y = c_vector_alloc (a->size1);
	n = (long) a->size1;
	m = (long) a->size2;
	lda = (long) a->lda;
	incx = (long) x->stride;
	incy = (long) y->stride;
	dgemv_ (&trans, &n, &m, &alpha, a->data, &lda, x->data, &incx, &beta, y->data, &incy);
	return y;
}

/* y = alpha * a' * x + beta */
c_vector *
c_matrix_transpose_dot_vector (double alpha, const c_matrix *a, const c_vector *x, double beta)
{
	char		trans = 'T';
	long		n;
	long		m;
	long		lda;
	long		incx;
	long		incy;
	c_vector	*y;
	if (c_matrix_is_empty (a)) c_error ("c_matrix_transpose_dot_vector", "matrix is empty.");
	if (c_vector_is_empty (x)) c_error ("c_matrix_transpose_dot_vector", "vector is empty.");
	if (a->size1 != x->size) c_error ("c_matrix_transpose_dot_vector", "vector and matrix size does not match.");
	y = c_vector_alloc (a->size2);
	n = (long) a->size1;
	m = (long) a->size2;
	lda = (long) a->lda;
	incx = (long) x->stride;
	incy = (long) y->stride;
	dgemv_ (&trans, &n, &m, &alpha, a->data, &lda, x->data, &incx, &beta, y->data, &incy);
	return y;
}

/* c = alpha * a * b + beta */
c_matrix *
c_matrix_dot_matrix (double alpha, const c_matrix *a, const c_matrix *b, double beta)
{
	char		transA = 'N';
	char		transB = 'N';
	long		m;
	long		n;
	long		k;
	long		lda;
	long		ldb;
	long		ldc;
	c_matrix	*c;
	if (c_matrix_is_empty (a)) c_error ("c_matrix_dot_matrix", "matrix *a is empty.");
	if (c_matrix_is_empty (b)) c_error ("c_matrix_dot_matrix", "matrix *b is empty.");
	if (a->size1 != b->size1) c_error ("c_matrix_dot_matrix", "matrix size does not match.");
	c = c_matrix_alloc (a->size1, b->size2);
	m = (long) a->size1;
	n = (long) b->size2;
	k = (long) a->size2;
	lda = (long) a->lda;
	ldb = (long) b->lda;
	ldc = (long) c->lda;
	dgemm_ (&transA, &transB, &m, &n, &k, &alpha, a->data, &lda, b->data, &ldb, &beta, c->data, &ldc);
	return c;
}

/* c = alpha * a * b' + beta */
c_matrix *
c_matrix_dot_matrix_transpose (double alpha, const c_matrix *a, const c_matrix *b, double beta)
{
	char		transA = 'N';
	char		transB = 'T';
	long		m;
	long		n;
	long		k;
	long		lda;
	long		ldb;
	long		ldc;
	c_matrix	*c;
	if (c_matrix_is_empty (a)) c_error ("c_matrix_dot_matrix_transpose", "matrix *a is empty.");
	if (c_matrix_is_empty (b)) c_error ("c_matrix_dot_matrix_transpose", "matrix *b is empty.");
	if (a->size1 != b->size1) c_error ("c_matrix_dot_matrix_transpose", "matrix size does not match.");
	c = c_matrix_alloc (a->size2, b->size1);
	m = (long) a->size1;
	n = (long) b->size1;
	k = (long) a->size2;
	lda = (long) a->lda;
	ldb = (long) b->lda;
	ldc = (long) c->lda;
	dgemm_ (&transA, &transB, &m, &n, &k, &alpha, a->data, &lda, b->data, &ldb, &beta, c->data, &ldc);
	return c;
}

/* c = alpha * a' * b + beta */
c_matrix *
c_matrix_transpose_dot_matrix (double alpha, const c_matrix *a, const c_matrix *b, double beta)
{
	char		transA = 'T';
	char		transB = 'N';
	long		m;
	long		n;
	long		k;
	long		lda;
	long		ldb;
	long		ldc;
	c_matrix	*c;
	if (c_matrix_is_empty (a)) c_error ("c_matrix_transpose_dot_matrix", "matrix *a is empty.");
	if (c_matrix_is_empty (b)) c_error ("c_matrix_transpose_dot_matrix", "matrix *b is empty.");
	if (a->size1 != b->size1) c_error ("c_matrix_transpose_dot_matrix", "matrix size does not match.");
	c = c_matrix_alloc (a->size2, b->size2);
	m = (long) a->size2;
	n = (long) b->size2;
	k = (long) a->size1;
	lda = (long) a->lda;
	ldb = (long) b->lda;
	ldc = (long) c->lda;
	dgemm_ (&transA, &transB, &m, &n, &k, &alpha, a->data, &lda, b->data, &ldb, &beta, c->data, &ldc);
	return c;
}
