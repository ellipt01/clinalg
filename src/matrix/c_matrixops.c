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
extern void	daxpy_ (long *n, double *alpha, double *x, long *incx, double *y, long *incy);
extern void	dgemv_ (char *trans, long *n, long *m, double *alpha, double *a, long *lda, double *x, long *incx, double *beta, double *y, long *incy);
extern void	dgemm_ (char *transA, char *transB, long *m, long *n, long *k, double *alpha, double *a, long *lda, double *b, long *ldb, double *beta, double *c, long *ldc);

/* lapack */
extern double	dlange_ (char *norm, long *m, long *n, double *data, long *lda, double *w);

/* x = x - y */
void
c_matrix_sub (c_matrix *x, const c_matrix *y)
{
	long	n;
	long	incx = 1;
	long	incy = 1;
	double	alpha = -1.;
	if (x->size1 != y->size1 || x->size2 != y->size2) c_error ("c_matrix_sub", "matrix size done not match.");

	if (x->size1 == x->lda || y->size1 == y->lda) {
		n = x->tsize;
		daxpy_ (&n, &alpha, y->data, &incx, x->data, &incy);
	} else {
		int		j;
		n = x->size1;
		for (j = 0; j < x->size2; j++) {
			double	*xj = x->data + INDEX_OF_MATRIX (x, 0, j);
			double	*yj = y->data + INDEX_OF_MATRIX (y, 0, j);
			daxpy_ (&n, &alpha, yj, &incx, xj, &incy);
		}
	}
	return;
}

double
c_matrix_nrm (c_matrix *a, char norm)
{
	long	m, n, lda;
	double	val;
	double	*w = NULL;

	m = a->size1;
	n = a->size2;
	lda = a->lda;

	switch (norm) {
		/* norm2 (A) */
		case '2':
//		val = _matrix_nrm2 (a);
//		return val;

		/* max (abs (A(i, j)) ) */
		case 'M':
		case 'm':

		/* norm1 (A) */
		case '1':
		case 'O':
		case 'o':

		/* normF (A) */
		case 'F':
		case 'f':
		case 'E':
		case 'e':
		break;

		/* normI (A) */
		case 'I':
		case 'i':
		w = (double *) malloc (m * sizeof (double));
		break;

		default:
		c_error ("c_matrix_nrm", "invalid norm.");
		break;
	}

	val = dlange_ (&norm, &m, &n, a->data, &lda, w);
	if (w) free (w);

	return val;
}

void
c_matrix_upper_triangular_memcpy (c_matrix *tr, const c_matrix *a)
{
	int			j;
	long		incx = 1;
	long		incy = 1;

	size_t		min_m = C_MIN (tr->size1, a->size1);
	size_t		min_n = C_MIN (tr->size2, a->size2);
	for (j = 0; j < min_n; j++) {
		long	n = (j + 1 < min_m) ? (long) (j + 1) : (long) min_m;
		dcopy_ (&n, a->data + j * a->lda, &incx, tr->data + j * tr->lda, &incy);
	}
	return;
}

void
c_matrix_lower_triangular_memcpy (c_matrix *tr, const c_matrix *a)
{
	int			j;
	long		incx = 1;
	long		incy = 1;

	size_t		min_m = C_MIN (tr->size1, a->size1);
	size_t		min_n = C_MIN (tr->size2, a->size2);
	for (j = 0; j < min_n; j++) {
		long	n;
		n = min_m - j;
		if (n <= 0) break;
		dcopy_ (&n, a->data + j * (a->lda + 1), &incx, tr->data + j * (tr->lda + 1), &incy);
	}
	return;
}

c_vector *
c_matrix_get_diagonal (const c_matrix *a)
{
	long		n;
	long		lda;
	long		stride;
	size_t		min_mn = C_MIN (a->size1, a->size2);
	c_vector	*d = c_vector_alloc (min_mn);

	n = (long) min_mn;
	lda = (long) a->lda;
	stride = (long) d->stride;
	dcopy_ (&n, a->data, &lda, d->data, &stride);

	return d;
}

void
c_matrix_set_diagonal (const c_vector *d, c_matrix *a)
{
	long		n;
	long		stride;
	long		lda;
	size_t		min_mn = C_MIN (a->size1, a->size2);

	n = (long) C_MIN (d->size, min_mn);
	stride = d->stride;
	lda = a->lda + 1;
	dcopy_ (&n, d->data, &stride, a->data, &lda);

	return;
}

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
