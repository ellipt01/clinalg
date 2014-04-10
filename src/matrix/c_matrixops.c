/*
 * c_matrixops.c
 *
 *  Created on: 2014/04/03
 *      Author: utsugi
 */

#include <c_matrix.h>

extern void	c_error (const char * function_name, const char *error_msg);

/* blas */
extern void	dgemv_ (char *trans, long *n, long *m, double *alpha, double *a, long *lda, double *x, long *incx, double *beta, double *y, long *incy);

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

