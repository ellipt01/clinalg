/*
 * c_matrixops.c
 *
 *  Created on: 2014/04/03
 *      Author: utsugi
 */

#include <clinalg_macros.h>
#include <c_matrix.h>

/* c_linalg_util.c */
extern void	c_error (const char * function_name, const char *error_msg);

/* c_linalg_sv.c */
extern int		c_linalg_lapack_dgesvd (char jobu, char jobvt, c_matrix *a, c_matrix **u, c_matrix **vt, c_vector **s);

/* blas */
#ifndef HAVE_BLAS_H
extern void	dcopy_ (int *n, double *x, int *incx, double *y, int *incy);
extern void	daxpy_ (int *n, double *alpha, double *x, int *incx, double *y, int *incy);
extern void	dgemv_ (char *trans, int *n, int *m, double *alpha, double *a, int *lda, double *x, int *incx, double *beta, double *y, int *incy);
extern void	dsymv_ (char *uplo, int *n, double *alpha, double *a, int *lda, double *x, int *incx, double *beta, double *y, int *incy);
extern void	dgemm_ (char *transA, char *transB, int *m, int *n, int *k, double *alpha, double *a, int *lda, double *b, int *ldb, double *beta, double *c, int *ldc);
extern void	dsymm_ (char *side, char *uplo, int *m, int *n, double *alpha, double *a, int *lda, double *b, int *ldb, double *beta, double *c, int *ldc);
#endif

/* lapack */
#ifndef HAVE_LAPACK_H
extern double	dlange_ (char *norm, int *m, int *n, double *data, int *lda, double *w);
extern void	dswap_ (int *n, double *x, int *incx, double *y, int *incy);
#endif

/* y = x + y */
void
c_matrix_add (c_matrix *y, const c_matrix *x)
{
	int		n;
	int		incx = 1;
	int		incy = 1;
	double	alpha = 1.;
	if (x->size1 != y->size1 || x->size2 != y->size2) c_error ("c_matrix_add", "matrix size done not match.");

	if (x->size1 == x->lda || y->size1 == y->lda) {
		n = x->tsize;
		daxpy_ (&n, &alpha, x->data, &incx, y->data, &incy);
	} else {
		int		j;
		n = x->size1;
		for (j = 0; j < x->size2; j++) {
			double	*xj = POINTER_OF_MATRIX (x, 0, j);
			double	*yj = POINTER_OF_MATRIX (y, 0, j);
			daxpy_ (&n, &alpha, xj, &incx, yj, &incy);
		}
	}
	return;
}

/* y = y - x */
void
c_matrix_sub (c_matrix *y, const c_matrix *x)
{
	int		n;
	int		incx = 1;
	int		incy = 1;
	double	alpha = -1.;
	if (x->size1 != y->size1 || x->size2 != y->size2) c_error ("c_matrix_sub", "matrix size done not match.");

	if (x->size1 == x->lda || y->size1 == y->lda) {
		n = x->tsize;
		daxpy_ (&n, &alpha, x->data, &incx, y->data, &incy);
	} else {
		int		j;
		n = x->size1;
		for (j = 0; j < x->size2; j++) {
			double	*xj = POINTER_OF_MATRIX (x, 0, j);
			double	*yj = POINTER_OF_MATRIX (y, 0, j);
			daxpy_ (&n, &alpha, xj, &incx, yj, &incy);
		}
	}
	return;
}

/* x = x - y */
void
c_matrix_axpy (double alpha, const c_matrix *x, c_matrix *y)
{
	int		n;
	int		incx = 1;
	int		incy = 1;
	if (x->size1 != y->size1 || x->size2 != y->size2) c_error ("c_matrix_axpy", "matrix size done not match.");

	if (x->size1 == x->lda || y->size1 == y->lda) {
		n = x->tsize;
		daxpy_ (&n, &alpha, x->data, &incx, y->data, &incy);
	} else {
		int		j;
		n = x->size1;
		for (j = 0; j < x->size2; j++) {
			double	*xj = POINTER_OF_MATRIX (x, 0, j);
			double	*yj = POINTER_OF_MATRIX (y, 0, j);
			daxpy_ (&n, &alpha, xj, &incx, yj, &incy);
		}
	}
	return;
}

static double
c_matrix_nrm2 (const c_matrix *a)
{
	int			info;
	double		nrm = 0.;
	c_vector	*s;
	c_matrix	*tmp = c_matrix_alloc (a->size1, a->size2);
	c_matrix_memcpy (tmp, a);

	info = c_linalg_lapack_dgesvd ('N', 'N', tmp, NULL, NULL, &s);
	c_matrix_free (tmp);

	if (info == 0) nrm = c_vector_get (s, 0);
	c_vector_free (s);

	return nrm;
}

double
c_matrix_nrm (c_matrix *a, char norm)
{
	int		m, n, lda;
	double	val;
	double	*w = NULL;

	if (c_matrix_is_empty (a)) c_error ("c_matrix_nrm", "matrix is empty.");

	m = a->size1;
	n = a->size2;
	lda = a->lda;

	switch (norm) {
		/* norm2 (A) */
		case '2':
			val = c_matrix_nrm2 (a);
			return val;

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
c_matrix_swap_rows (const size_t i, const size_t j, c_matrix *a)
{
	int			n;
	int			inc;
	double		*rowi;
	double		*rowj;

	if (c_matrix_is_empty (a)) c_error ("c_matrix_swap_rows", "matrix is empty.");
	if (i < 0 || a->size1 <= i) c_error ("c_matrix_swap_rows", "first index out of range.");
	if (j < 0 || a->size1 <= j) c_error ("c_matrix_swap_rows", "second index out of range.");
	if (i == j) return;

	n = (int) a->size2;
	inc = (int) a->lda;
	rowi = POINTER_OF_MATRIX (a, i, 0);
	rowj = POINTER_OF_MATRIX (a, j, 0);

	dswap_ (&n, rowi, &inc, rowj, &inc);

	return;
}

void
c_matrix_swap_cols (const size_t i, const size_t j, c_matrix *a)
{
	int			n;
	int			inc;
	double		*coli;
	double		*colj;

	if (c_matrix_is_empty (a)) c_error ("c_matrix_swap_rows", "matrix is empty.");
	if (i < 0 || a->size2 <= i) c_error ("c_matrix_swap_rows", "first index out of range.");
	if (j < 0 || a->size2 <= j) c_error ("c_matrix_swap_rows", "second index out of range.");
	if (i == j) return;

	n = (int) a->size1;
	inc = 1;
	coli = POINTER_OF_MATRIX (a, 0, i);
	colj = POINTER_OF_MATRIX (a, 0, j);

	dswap_ (&n, coli, &inc, colj, &inc);

	return;
}

void
c_matrix_permute_rows (c_matrix *a, const c_vector_int *p)
{
	int			i;
	size_t		size;

	if (c_matrix_is_empty (a)) c_error ("c_matrix_permute_rows", "matrix is empty.");
	if (c_vector_int_is_empty (p)) c_error ("c_matrix_permute_rows", "permutation is empty.");

	size = (size_t) C_MIN (a->size1, p->size);
	for (i = 0; i < size; i++) {
		int		pi = c_vector_int_get (p, i) - 1;
		if (pi < 0 || pi == i || a->size1 <= pi) continue;
		c_matrix_swap_rows (i, pi, a);
	}
	return;
}

void
c_matrix_permute_cols (c_matrix *a, const c_vector_int *p)
{
	int			i;
	size_t		size;

	if (c_matrix_is_empty (a)) c_error ("c_matrix_permute_cols", "matrix is empty.");
	if (c_vector_int_is_empty (p)) c_error ("c_matrix_permute_cols", "permutation is empty.");

	size = (size_t) C_MIN (a->size2, p->size);
	for (i = 0; i < size; i++) {
		int		pi = c_vector_int_get (p, i) - 1;
		if (pi < 0 || pi == i || a->size2 <= pi) continue;
		c_matrix_swap_cols (i, pi, a);
	}
	return;
}

void
c_matrix_add_rowcols (c_matrix *a, const size_t dm, const size_t dn)
{
	int			i, j;
	int			n;
	int			inc = 1;
	size_t		lda;

	if (c_matrix_is_empty (a)) c_error ("c_matrix_add_rowcols", "matrix is empty.");
	if (!a->owner) c_error ("c_matrix_add_rowcols", "cannot resize matrix_view.");

	if (dm <= 0 && dn <= 0) return;

	n = (int) a->size1;
	lda = (int) a->lda;
	if (dm > 0) {
		a->size1 += dm;
		a->lda += dm;
	}
	if (dn > 0) a->size2 += dn;
	a->tsize = a->lda * a->size2;
	a->data = (double *) realloc (a->data, a->tsize * sizeof (double));

	if (dm > 0) {
		c_vector	*col = c_vector_alloc (n);
		for (j = (a->size2 - dn) - 1; 0 < j; j--) {
			dcopy_ (&n, a->data + j * lda, &inc, col->data, &inc);
			dcopy_ (&n, col->data, &inc, POINTER_OF_MATRIX (a, 0, j), &inc);
		}
		c_vector_free (col);
	}
	for (i = a->size1 - dm; i < a->size1; i++) {
		for (j = 0; j < a->size2; j++) c_matrix_set (a, i, j, 0.);
	}
	for (j = a->size2 - dn; j < a->size2; j++) {
		for (i = 0; i < a->size1 - dm; i++) c_matrix_set (a, i, j, 0.);
	}

	return;
}

void
c_matrix_remove_rowcols (c_matrix *a, const size_t dm, const size_t dn)
{
	int			j;
	int			n;
	int			inc = 1;
	size_t		lda;

	if (c_matrix_is_empty (a)) c_error ("c_matrix_remove_rowcols", "matrix is empty.");
	if (!a->owner) c_error ("c_matrix_remove_rowcols", "cannot resize matrix_view.");
	if (a->size1 - dm < 0) c_error ("c_matrix_remove_rowcols", "a->size1 must be > dm.");
	if (a->size2 - dn < 0) c_error ("c_matrix_remove_rowcols", "a->size2 must be > dn.");

	if (dm <= 0 && dn <= 0) return;

	lda = (int) a->lda;
	if (dm > 0) {
		a->size1 -= dm;
		a->lda -= dm;
	}
	if (dn > 0) a->size2 -= dn;
	a->tsize = a->lda * a->size2;
	n = (int) a->size1;

	if (dm > 0 && a->size1 > 0 && a->size2 > 0) {
		c_vector	*col = c_vector_alloc (a->size1);
		for (j = 1; j < a->size2; j++) {
			dcopy_ (&n, a->data + j * lda, &inc, col->data, &inc);
			dcopy_ (&n, col->data, &inc, POINTER_OF_MATRIX (a, 0, j), &inc);
		}
		c_vector_free (col);
	}
	a->data = (double *) realloc (a->data, a->tsize * sizeof (double));

	return;
}

void
c_matrix_merge_row (c_matrix *a, const c_matrix *b)
{
	int			j, k;
	int			n;
	int			inc;
	if (c_matrix_is_empty (a)) c_error ("c_matrix_merge_row", "matrix *a is empty.");
	if (c_matrix_is_empty (b)) c_error ("c_matrix_merge_row", "matrix *b is empty.");
	if (a->size1 != b->size1) c_error ("c_matrix_merge_row", "matrix size dose not match.");

	k = (int) a->size2;
	c_matrix_add_rowcols (a, 0, b->size2);

	n = (int) b->size1;
	inc = 1;
	for (j = 0; j < b->size2; j++) {
		double	*bj = POINTER_OF_MATRIX (b, 0, j);
		dcopy_ (&n, bj, &inc, POINTER_OF_MATRIX (a, 0, k + j), &inc);
	}
	return;
}

void
c_matrix_merge_col (c_matrix *a, const c_matrix *b)
{
	int			j, k;
	int			n;
	int			inc;
	if (c_matrix_is_empty (a)) c_error ("c_matrix_merge_col", "matrix *a is empty.");
	if (c_matrix_is_empty (b)) c_error ("c_matrix_merge_col", "matrix *b is empty.");
	if (a->size2 != b->size2) c_error ("c_matrix_merge_col", "matrix size dose not match.");

	k = (int) a->size1;
	c_matrix_add_rowcols (a, b->size1, 0);

	n = (int) b->size1;
	inc = 1;
	for (j = 0; j < b->size2; j++) {
		double	*bj = POINTER_OF_MATRIX (b, 0, j);
		dcopy_ (&n, bj, &inc, POINTER_OF_MATRIX (a, k, j), &inc);
	}
	return;
}

c_matrix *
c_matrix_identity (const size_t size1, const size_t size2)
{
	int			i;
	size_t		min_mn = C_MIN (size1, size2);
	c_matrix	*c = c_matrix_alloc (size1, size2);
	c_matrix_set_zero (c);
	for (i = 0; i < min_mn; i++) c->data[i * (c->lda + 1)] = 1.;
	return c;
}

c_matrix *
c_matrix_transpose (c_matrix *a)
{
	int			j;
	c_vector	*col;
	c_matrix	*at = c_matrix_alloc (a->size2, a->size1);

	int			n;
	int			incx;
	int			incy = (int) at->lda;
	for (j = 0; j < a->size2; j++) {
		col = c_matrix_column (a, j);
		n = (int) a->size1;
		incx = (int) col->stride;
		dcopy_ (&n, col->data, &incx, at->data + j, &incy);
	}
	return at;
}

/* y = alpha * a * x + beta */
c_vector *
c_matrix_dot_vector (double alpha, const c_matrix *a, const c_vector *x, double beta)
{
	char		trans = 'N';
	int			n;
	int			m;
	int			lda;
	int			incx;
	int			incy;
	c_vector	*y;

	if (c_matrix_is_empty (a)) c_error ("c_matrix_dot_vector", "matrix is empty.");
	if (c_vector_is_empty (x)) c_error ("c_matrix_dot_vector", "vector is empty.");
	if (a->size2 != x->size) c_error ("c_matrix_dot_vector", "vector and matrix size does not match.");

	y = c_vector_alloc (a->size1);
	n = (int) a->size1;
	m = (int) a->size2;
	lda = (int) a->lda;
	incx = (int) x->stride;
	incy = (int) y->stride;
	dgemv_ (&trans, &n, &m, &alpha, a->data, &lda, x->data, &incx, &beta, y->data, &incy);
	return y;
}

/* y = alpha * a' * x + beta */
c_vector *
c_matrix_transpose_dot_vector (double alpha, const c_matrix *a, const c_vector *x, double beta)
{
	char		trans = 'T';
	int			n;
	int			m;
	int			lda;
	int			incx;
	int			incy;
	c_vector	*y;

	if (c_matrix_is_empty (a)) c_error ("c_matrix_transpose_dot_vector", "matrix is empty.");
	if (c_vector_is_empty (x)) c_error ("c_matrix_transpose_dot_vector", "vector is empty.");
	if (a->size1 != x->size) c_error ("c_matrix_transpose_dot_vector", "vector and matrix size does not match.");

	y = c_vector_alloc (a->size2);
	n = (int) a->size1;
	m = (int) a->size2;
	lda = (int) a->lda;
	incx = (int) x->stride;
	incy = (int) y->stride;
	dgemv_ (&trans, &n, &m, &alpha, a->data, &lda, x->data, &incx, &beta, y->data, &incy);
	return y;
}

/* y = alpha * a * x + beta, only refer upper triangular part of a */
c_vector *
c_matrix_symm_upper_dot_vector (double alpha, const c_matrix *a, const c_vector *x, double beta)
{
	int			n;
	int			lda;
	int			incx;
	int			incy;
	c_vector	*y;

	if (c_matrix_is_empty (a)) c_error ("c_matrix_symm_upper_dot_vector", "matrix is empty.");
	if (c_vector_is_empty (x)) c_error ("c_matrix_symm_upper_dot_vector", "vector is empty.");
	if (!c_matrix_is_square (a)) c_error ("c_matrix_symm_upper_dot_vector", "vector and matrix size does not match.");

	y = c_vector_alloc (a->size1);
	n = (int) a->size1;
	lda = (int) a->lda;
	incx = (int) x->stride;
	incy = (int) y->stride;
	dsymv_ ("U", &n, &alpha, a->data, &lda, x->data, &incx, &beta, y->data, &incy);
	return y;
}

/* y = alpha * a * x + beta, only refer lower triangular part of a */
c_vector *
c_matrix_symm_lower_dot_vector (double alpha, const c_matrix *a, const c_vector *x, double beta)
{
	int			n;
	int			lda;
	int			incx;
	int			incy;
	c_vector	*y;

	if (c_matrix_is_empty (a)) c_error ("c_matrix_symm_lower_dot_vector", "matrix is empty.");
	if (c_vector_is_empty (x)) c_error ("c_matrix_symm_lower_dot_vector", "vector is empty.");
	if (!c_matrix_is_square (a)) c_error ("c_matrix_symm_lower_dot_vector", "vector and matrix size does not match.");

	y = c_vector_alloc (a->size1);
	n = (int) a->size1;
	lda = (int) a->lda;
	incx = (int) x->stride;
	incy = (int) y->stride;
	dsymv_ ("L", &n, &alpha, a->data, &lda, x->data, &incx, &beta, y->data, &incy);
	return y;
}

/* c = alpha * a * b + beta */
c_matrix *
c_matrix_dot_matrix (double alpha, const c_matrix *a, const c_matrix *b, double beta)
{
	char		transA = 'N';
	char		transB = 'N';
	int			m;
	int			n;
	int			k;
	int			lda;
	int			ldb;
	int			ldc;
	c_matrix	*c;

	if (c_matrix_is_empty (a)) c_error ("c_matrix_dot_matrix", "matrix *a is empty.");
	if (c_matrix_is_empty (b)) c_error ("c_matrix_dot_matrix", "matrix *b is empty.");
	if (a->size2 != b->size1) c_error ("c_matrix_dot_matrix", "matrix size does not match.");

	c = c_matrix_alloc (a->size1, b->size2);
	m = (int) a->size1;
	n = (int) b->size2;
	k = (int) a->size2;
	lda = (int) a->lda;
	ldb = (int) b->lda;
	ldc = (int) c->lda;
	dgemm_ (&transA, &transB, &m, &n, &k, &alpha, a->data, &lda, b->data, &ldb, &beta, c->data, &ldc);
	return c;
}

/* c = alpha * a * b' + beta */
c_matrix *
c_matrix_dot_matrix_transpose (double alpha, const c_matrix *a, const c_matrix *b, double beta)
{
	char		transA = 'N';
	char		transB = 'T';
	int			m;
	int			n;
	int			k;
	int			lda;
	int			ldb;
	int			ldc;
	c_matrix	*c;

	if (c_matrix_is_empty (a)) c_error ("c_matrix_dot_matrix_transpose", "matrix *a is empty.");
	if (c_matrix_is_empty (b)) c_error ("c_matrix_dot_matrix_transpose", "matrix *b is empty.");
	if (a->size2 != b->size2) c_error ("c_matrix_dot_matrix_transpose", "matrix size does not match.");

	c = c_matrix_alloc (a->size1, b->size1);
	m = (int) a->size1;
	n = (int) b->size1;
	k = (int) a->size2;
	lda = (int) a->lda;
	ldb = (int) b->lda;
	ldc = (int) c->lda;
	dgemm_ (&transA, &transB, &m, &n, &k, &alpha, a->data, &lda, b->data, &ldb, &beta, c->data, &ldc);
	return c;
}

/* c = alpha * a' * b + beta */
c_matrix *
c_matrix_transpose_dot_matrix (double alpha, const c_matrix *a, const c_matrix *b, double beta)
{
	char		transA = 'T';
	char		transB = 'N';
	int			m;
	int			n;
	int			k;
	int			lda;
	int			ldb;
	int			ldc;
	c_matrix	*c;

	if (c_matrix_is_empty (a)) c_error ("c_matrix_transpose_dot_matrix", "matrix *a is empty.");
	if (c_matrix_is_empty (b)) c_error ("c_matrix_transpose_dot_matrix", "matrix *b is empty.");
	if (a->size1 != b->size1) c_error ("c_matrix_transpose_dot_matrix", "matrix size does not match.");

	c = c_matrix_alloc (a->size2, b->size2);
	m = (int) a->size2;
	n = (int) b->size2;
	k = (int) a->size1;
	lda = (int) a->lda;
	ldb = (int) b->lda;
	ldc = (int) c->lda;
	dgemm_ (&transA, &transB, &m, &n, &k, &alpha, a->data, &lda, b->data, &ldb, &beta, c->data, &ldc);
	return c;
}

/* c = alpha * a' * b' + beta */
c_matrix *
c_matrix_transpose_dot_matrix_transpose (double alpha, const c_matrix *a, const c_matrix *b, double beta)
{
	char		transA = 'T';
	char		transB = 'T';
	int			m;
	int			n;
	int			k;
	int			lda;
	int			ldb;
	int			ldc;
	c_matrix	*c;

	if (c_matrix_is_empty (a)) c_error ("c_matrix_transpose_dot_matrix_transpose", "matrix *a is empty.");
	if (c_matrix_is_empty (b)) c_error ("c_matrix_transpose_dot_matrix_transpose", "matrix *b is empty.");
	if (a->size1 != b->size2) c_error ("c_matrix_transpose_dot_matrix_transpose", "matrix size does not match.");

	c = c_matrix_alloc (a->size2, b->size1);
	m = (int) a->size2;
	n = (int) b->size1;
	k = (int) a->size1;
	lda = (int) a->lda;
	ldb = (int) b->lda;
	ldc = (int) c->lda;
	dgemm_ (&transA, &transB, &m, &n, &k, &alpha, a->data, &lda, b->data, &ldb, &beta, c->data, &ldc);
	return c;
}

/* c = alpha * a * b + beta */
static c_matrix *
c_matrix_dot_matrix_symm (char *side, char *uplo, double alpha, const c_matrix *a, const c_matrix *b, double beta)
{
	int			m;
	int			n;
	int			k;
	int			lda;
	int			ldb;
	int			ldc;
	c_matrix	*c = c_matrix_alloc (a->size2, b->size1);
	m = (int) a->size2;
	n = (int) b->size1;
	k = (int) a->size1;
	lda = (int) a->lda;
	ldb = (int) b->lda;
	ldc = (int) c->lda;
	dsymm_ (side, uplo, &m, &n, &alpha, a->data, &lda, b->data, &ldb, &beta, c->data, &ldc);
	return c;
}

c_matrix *
c_matrix_symm_upper_dot_matrix (double alpha, const c_matrix *a, const c_matrix *b, double beta)
{
	if (c_matrix_is_empty (a)) c_error ("c_matrix_symm_upper_dot_matrix", "matrix *a is empty.");
	if (c_matrix_is_empty (b)) c_error ("c_matrix_symm_upper_dot_matrix", "matrix *b is empty.");
	if (a->size1 != b->size2) c_error ("c_matrix_symm_upper_dot_matrix", "matrix size does not match.");

	return c_matrix_dot_matrix_symm ("L", "U", alpha, a, b, beta);
}

c_matrix *
c_matrix_symm_lower_dot_matrix (double alpha, const c_matrix *a, const c_matrix *b, double beta)
{
	if (c_matrix_is_empty (a)) c_error ("c_matrix_symm_lower_dot_matrix", "matrix *a is empty.");
	if (c_matrix_is_empty (b)) c_error ("c_matrix_symm_lower_dot_matrix", "matrix *b is empty.");
	if (a->size1 != b->size2) c_error ("c_matrix_symm_lower_dot_matrix", "matrix size does not match.");

	return c_matrix_dot_matrix_symm ("L", "L", alpha, a, b, beta);
}

c_matrix *
c_matrix_dot_matrix_symm_upper (double alpha, const c_matrix *a, const c_matrix *b, double beta)
{
	if (c_matrix_is_empty (a)) c_error ("c_matrix_dot_matrix_symm_upper", "matrix *a is empty.");
	if (c_matrix_is_empty (b)) c_error ("c_matrix_dot_matrix_symm_upper", "matrix *b is empty.");
	if (a->size1 != b->size2) c_error ("c_matrix_dot_matrix_symm_upper", "matrix size does not match.");

	return c_matrix_dot_matrix_symm ("R", "U", alpha, a, b, beta);
}

c_matrix *
c_matrix_dot_matrix_symm_lower (double alpha, const c_matrix *a, const c_matrix *b, double beta)
{
	if (c_matrix_is_empty (a)) c_error ("c_matrix_dot_matrix_symm_lower", "matrix *a is empty.");
	if (c_matrix_is_empty (b)) c_error ("c_matrix_symm_lower_dot_matrix", "matrix *b is empty.");
	if (a->size1 != b->size2) c_error ("c_matrix_dot_matrix_symm_lower", "matrix size does not match.");

	return c_matrix_dot_matrix_symm ("R", "L", alpha, a, b, beta);
}
