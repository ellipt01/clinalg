/*
 * c_matrixops.c
 *
 *  Created on: 2014/04/03
 *      Author: utsugi
 */

#include <c_matrix.h>
#include <c_linalg_lapack.h>
#include "../../include/clinalg_macros.h"
#include "../../include/clinalg_utils.h"

#include "private.h"

/* y = x + y */
void
c_matrix_add (c_matrix *y, const c_matrix *x)
{
	double	alpha = 1.;
	if (x->size1 != y->size1 || x->size2 != y->size2) c_error ("c_matrix_add", "matrix size done not match.");

	if (x->size1 == x->lda || y->size1 == y->lda) F77CALL (daxpy) (&x->tsize, &alpha, x->data, &ione, y->data, &ione);
	else {
		int		j;
		for (j = 0; j < x->size2; j++) {
			double	*xj = POINTER_OF_MATRIX (x, 0, j);
			double	*yj = POINTER_OF_MATRIX (y, 0, j);
			F77CALL (daxpy) (&x->size1, &alpha, xj, &ione, yj, &ione);
		}
	}
	return;
}

/* y = y - x */
void
c_matrix_sub (c_matrix *y, const c_matrix *x)
{
	double	alpha = -1.;
	if (x->size1 != y->size1 || x->size2 != y->size2) c_error ("c_matrix_sub", "matrix size done not match.");

	if (x->size1 == x->lda || y->size1 == y->lda) F77CALL (daxpy) (&x->tsize, &alpha, x->data, &ione, y->data, &ione);
	else {
		int		j;
		for (j = 0; j < x->size2; j++) {
			double	*xj = POINTER_OF_MATRIX (x, 0, j);
			double	*yj = POINTER_OF_MATRIX (y, 0, j);
			F77CALL (daxpy) (&x->size1, &alpha, xj, &ione, yj, &ione);
		}
	}
	return;
}

/* x = x - y */
void
c_matrix_axpy (const double alpha, const c_matrix *x, c_matrix *y)
{
	if (x->size1 != y->size1 || x->size2 != y->size2) c_error ("c_matrix_axpy", "matrix size done not match.");

	if (x->size1 == x->lda || y->size1 == y->lda) F77CALL (daxpy) (&x->tsize, &alpha, x->data, &ione, y->data, &ione);
	else {
		int		j;
		for (j = 0; j < x->size2; j++) {
			double	*xj = POINTER_OF_MATRIX (x, 0, j);
			double	*yj = POINTER_OF_MATRIX (y, 0, j);
			F77CALL (daxpy) (&x->size1, &alpha, xj, &ione, yj, &ione);
		}
	}
	return;
}

/* x = alpha * x */
void
c_matrix_scale (const double alpha, const c_matrix *x)
{
	F77CALL (dscal) (&x->tsize, &alpha, x->data, &ione);
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
c_matrix_nrm (c_matrix *a, const char norm)
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

	val = F77CALL (dlange) (&norm, &m, &n, a->data, &lda, w);
	if (w) free (w);

	return val;
}

void
c_matrix_swap_rows (const int i, const int j, c_matrix *a)
{
	double		*rowi;
	double		*rowj;

	if (c_matrix_is_empty (a))  c_error ("c_matrix_swap_rows", "matrix is empty.");
	if (i < 0 || a->size1 <= i) c_error ("c_matrix_swap_rows", "first index out of range.");
	if (j < 0 || a->size1 <= j) c_error ("c_matrix_swap_rows", "second index out of range.");
	if (i == j) return;

	rowi = POINTER_OF_MATRIX (a, i, 0);
	rowj = POINTER_OF_MATRIX (a, j, 0);

	F77CALL (dswap) (&a->size2, rowi, &a->lda, rowj, &a->lda);

	return;
}

void
c_matrix_swap_cols (const int i, const int j, c_matrix *a)
{
	double		*coli;
	double		*colj;

	if (c_matrix_is_empty (a))  c_error ("c_matrix_swap_rows", "matrix is empty.");
	if (i < 0 || a->size2 <= i) c_error ("c_matrix_swap_rows", "first index out of range.");
	if (j < 0 || a->size2 <= j) c_error ("c_matrix_swap_rows", "second index out of range.");
	if (i == j) return;

	coli = POINTER_OF_MATRIX (a, 0, i);
	colj = POINTER_OF_MATRIX (a, 0, j);

	F77CALL (dswap) (&a->size1, coli, &ione, colj, &ione);

	return;
}

void
c_matrix_permute_rows (c_matrix *a, const c_vector_int *p)
{
	int		i;
	int		size;

	if (c_matrix_is_empty (a)) c_error ("c_matrix_permute_rows", "matrix is empty.");
	if (c_vector_int_is_empty (p)) c_error ("c_matrix_permute_rows", "permutation is empty.");

	size = (int) C_MIN (a->size1, p->size);
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
	int		i;
	int		size;

	if (c_matrix_is_empty (a)) c_error ("c_matrix_permute_cols", "matrix is empty.");
	if (c_vector_int_is_empty (p)) c_error ("c_matrix_permute_cols", "permutation is empty.");

	size = (int) C_MIN (a->size2, p->size);
	for (i = 0; i < size; i++) {
		int		pi = c_vector_int_get (p, i) - 1;
		if (pi < 0 || pi == i || a->size2 <= pi) continue;
		c_matrix_swap_cols (i, pi, a);
	}
	return;
}

void
c_matrix_add_rowcols (c_matrix *a, const int dm, const int dn)
{
	int		i, j;
	int		n;
	int		lda;

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
			F77CALL (dcopy) (&n, a->data + j * lda, &ione, col->data, &ione);
			F77CALL (dcopy) (&n, col->data, &ione, POINTER_OF_MATRIX (a, 0, j), &ione);
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
c_matrix_remove_rowcols (c_matrix *a, const int dm, const int dn)
{
	int		j;
	int		n;
	int		lda;

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
			F77CALL (dcopy) (&n, a->data + j * lda, &ione, col->data, &ione);
			F77CALL (dcopy) (&n, col->data, &ione, POINTER_OF_MATRIX (a, 0, j), &ione);
		}
		c_vector_free (col);
	}
	a->data = (double *) realloc (a->data, a->tsize * sizeof (double));

	return;
}

void
c_matrix_merge_row (c_matrix *a, const c_matrix *b)
{
	int		j, k;
	if (c_matrix_is_empty (a)) c_error ("c_matrix_merge_row", "matrix *a is empty.");
	if (c_matrix_is_empty (b)) c_error ("c_matrix_merge_row", "matrix *b is empty.");
	if (a->size1 != b->size1)  c_error ("c_matrix_merge_row", "matrix size dose not match.");

	k = a->size2;
	c_matrix_add_rowcols (a, 0, b->size2);

	for (j = 0; j < b->size2; j++) {
		double	*bj = POINTER_OF_MATRIX (b, 0, j);
		F77CALL (dcopy) (&b->size1, bj, &ione, POINTER_OF_MATRIX (a, 0, k + j), &ione);
	}
	return;
}

void
c_matrix_merge_col (c_matrix *a, const c_matrix *b)
{
	int		j, k;
	if (c_matrix_is_empty (a)) c_error ("c_matrix_merge_col", "matrix *a is empty.");
	if (c_matrix_is_empty (b)) c_error ("c_matrix_merge_col", "matrix *b is empty.");
	if (a->size2 != b->size2)  c_error ("c_matrix_merge_col", "matrix size dose not match.");

	k = a->size1;
	c_matrix_add_rowcols (a, b->size1, 0);

	for (j = 0; j < b->size2; j++) {
		double	*bj = POINTER_OF_MATRIX (b, 0, j);
		F77CALL (dcopy) (&b->size1, bj, &ione, POINTER_OF_MATRIX (a, k, j), &ione);
	}
	return;
}

c_matrix *
c_matrix_identity (const int size1, const int size2)
{
	int			i;
	int			min_mn = C_MIN (size1, size2);
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
	for (j = 0; j < a->size2; j++) {
		col = c_matrix_column (a, j);
		F77CALL (dcopy) (&a->size1, col->data, &col->stride, at->data + j, &at->lda);
	}
	return at;
}

/* y = alpha * a * x */
c_vector *
c_matrix_dot_vector (const double alpha, const c_matrix *a, const c_vector *x)
{
	char		trans = 'N';
	c_vector	*y;

	if (c_matrix_is_empty (a)) c_error ("c_matrix_dot_vector", "matrix is empty.");
	if (c_vector_is_empty (x)) c_error ("c_matrix_dot_vector", "vector is empty.");
	if (a->size2 != x->size)   c_error ("c_matrix_dot_vector", "vector and matrix size does not match.");

	y = c_vector_alloc (a->size1);
	F77CALL (dgemv) (&trans, &a->size1, &a->size2, &alpha, a->data, &a->lda, x->data, &x->stride, &dzero, y->data, &y->stride);
	return y;
}

/* y = alpha * a' * x */
c_vector *
c_matrix_transpose_dot_vector (const double alpha, const c_matrix *a, const c_vector *x)
{
	char		trans = 'T';
	c_vector	*y;

	if (c_matrix_is_empty (a)) c_error ("c_matrix_transpose_dot_vector", "matrix is empty.");
	if (c_vector_is_empty (x)) c_error ("c_matrix_transpose_dot_vector", "vector is empty.");
	if (a->size1 != x->size)   c_error ("c_matrix_transpose_dot_vector", "vector and matrix size does not match.");

	y = c_vector_alloc (a->size2);
	F77CALL (dgemv) (&trans, &a->size1, &a->size2, &alpha, a->data, &a->lda, x->data, &x->stride, &dzero, y->data, &y->stride);
	return y;
}

/* y = alpha * a * x, only refer upper triangular part of a */
c_vector *
c_matrix_symm_upper_dot_vector (const double alpha, const c_matrix *a, const c_vector *x)
{
	c_vector	*y;

	if (c_matrix_is_empty (a))   c_error ("c_matrix_symm_upper_dot_vector", "matrix is empty.");
	if (c_vector_is_empty (x))   c_error ("c_matrix_symm_upper_dot_vector", "vector is empty.");
	if (!c_matrix_is_square (a)) c_error ("c_matrix_symm_upper_dot_vector", "vector and matrix size does not match.");

	y = c_vector_alloc (a->size1);
	F77CALL (dsymv) ("U", &a->size1, &alpha, a->data, &a->lda, x->data, &x->stride, &dzero, y->data, &y->stride);
	return y;
}

/* y = alpha * a * x, only refer lower triangular part of a */
c_vector *
c_matrix_symm_lower_dot_vector (const double alpha, const c_matrix *a, const c_vector *x)
{
	c_vector	*y;

	if (c_matrix_is_empty (a))   c_error ("c_matrix_symm_lower_dot_vector", "matrix is empty.");
	if (c_vector_is_empty (x))   c_error ("c_matrix_symm_lower_dot_vector", "vector is empty.");
	if (!c_matrix_is_square (a)) c_error ("c_matrix_symm_lower_dot_vector", "vector and matrix size does not match.");

	y = c_vector_alloc (a->size1);
	F77CALL (dsymv) ("L", &a->size1, &alpha, a->data, &a->lda, x->data, &x->stride, &dzero, y->data, &y->stride);
	return y;
}

/* c = alpha * a * b */
c_matrix *
c_matrix_dot_matrix (const double alpha, const c_matrix *a, const c_matrix *b)
{
	char		transA = 'N';
	char		transB = 'N';
	c_matrix	*c;

	if (c_matrix_is_empty (a)) c_error ("c_matrix_dot_matrix", "matrix *a is empty.");
	if (c_matrix_is_empty (b)) c_error ("c_matrix_dot_matrix", "matrix *b is empty.");
	if (a->size2 != b->size1)  c_error ("c_matrix_dot_matrix", "matrix size does not match.");

	c = c_matrix_alloc (a->size1, b->size2);
	F77CALL (dgemm) (&transA, &transB, &a->size1, &b->size2, &a->size2, &alpha, a->data, &a->lda, b->data, &b->lda, &dzero, c->data, &c->lda);
	return c;
}

/* c = alpha * a * b' */
c_matrix *
c_matrix_dot_matrix_transpose (const double alpha, const c_matrix *a, const c_matrix *b)
{
	char		transA = 'N';
	char		transB = 'T';
	c_matrix	*c;

	if (c_matrix_is_empty (a)) c_error ("c_matrix_dot_matrix_transpose", "matrix *a is empty.");
	if (c_matrix_is_empty (b)) c_error ("c_matrix_dot_matrix_transpose", "matrix *b is empty.");
	if (a->size2 != b->size2)  c_error ("c_matrix_dot_matrix_transpose", "matrix size does not match.");

	c = c_matrix_alloc (a->size1, b->size1);
	F77CALL (dgemm) (&transA, &transB, &a->size1, &b->size1, &a->size2, &alpha, a->data, &a->lda, b->data, &b->lda, &dzero, c->data, &c->lda);
	return c;
}

/* c = alpha * a' * b + beta */
c_matrix *
c_matrix_transpose_dot_matrix (const double alpha, const c_matrix *a, const c_matrix *b)
{
	char		transA = 'T';
	char		transB = 'N';
	c_matrix	*c;

	if (c_matrix_is_empty (a)) c_error ("c_matrix_transpose_dot_matrix", "matrix *a is empty.");
	if (c_matrix_is_empty (b)) c_error ("c_matrix_transpose_dot_matrix", "matrix *b is empty.");
	if (a->size1 != b->size1)  c_error ("c_matrix_transpose_dot_matrix", "matrix size does not match.");

	c = c_matrix_alloc (a->size2, b->size2);
	F77CALL (dgemm) (&transA, &transB, &a->size2, &b->size2, &a->size1, &alpha, a->data, &a->lda, b->data, &b->lda, &dzero, c->data, &c->lda);
	return c;
}

/* c = alpha * a' * b' + beta */
c_matrix *
c_matrix_transpose_dot_matrix_transpose (const double alpha, const c_matrix *a, const c_matrix *b)
{
	char		transA = 'T';
	char		transB = 'T';
	c_matrix	*c;

	if (c_matrix_is_empty (a)) c_error ("c_matrix_transpose_dot_matrix_transpose", "matrix *a is empty.");
	if (c_matrix_is_empty (b)) c_error ("c_matrix_transpose_dot_matrix_transpose", "matrix *b is empty.");
	if (a->size1 != b->size2)  c_error ("c_matrix_transpose_dot_matrix_transpose", "matrix size does not match.");

	c = c_matrix_alloc (a->size2, b->size1);
	F77CALL (dgemm) (&transA, &transB, &a->size2, &b->size1, &a->size1, &alpha, a->data, &a->lda, b->data, &b->lda, &dzero, c->data, &c->lda);
	return c;
}

/* c = alpha * a * b + beta */
static c_matrix *
c_matrix_dsymm (char *side, char *uplo, const double alpha, const c_matrix *a, const c_matrix *b)
{
	c_matrix	*c = c_matrix_alloc (a->size1, b->size2);
	F77CALL (dsymm) (side, uplo, &c->size1, &c->size2, &alpha, a->data, &a->lda, b->data, &b->lda, &dzero, c->data, &c->lda);
	return c;
}

c_matrix *
c_matrix_symm_upper_dot_matrix (const double alpha, const c_matrix *a, const c_matrix *b)
{
	if (c_matrix_is_empty (a)) c_error ("c_matrix_symm_upper_dot_matrix", "matrix *a is empty.");
	if (c_matrix_is_empty (b)) c_error ("c_matrix_symm_upper_dot_matrix", "matrix *b is empty.");
	if (a->size2 != b->size1)  c_error ("c_matrix_symm_upper_dot_matrix", "matrix size does not match.");

	return c_matrix_dsymm ("L", "U", alpha, a, b);
}

c_matrix *
c_matrix_symm_lower_dot_matrix (const double alpha, const c_matrix *a, const c_matrix *b)
{
	if (c_matrix_is_empty (a)) c_error ("c_matrix_symm_lower_dot_matrix", "matrix *a is empty.");
	if (c_matrix_is_empty (b)) c_error ("c_matrix_symm_lower_dot_matrix", "matrix *b is empty.");
	if (a->size2 != b->size1)  c_error ("c_matrix_symm_lower_dot_matrix", "matrix size does not match.");

	return c_matrix_dsymm ("L", "L", alpha, a, b);
}

c_matrix *
c_matrix_dot_matrix_symm_upper (const double alpha, const c_matrix *b, const c_matrix *a)
{
	if (c_matrix_is_empty (b)) c_error ("c_matrix_dot_matrix_symm_upper", "matrix *b is empty.");
	if (c_matrix_is_empty (a)) c_error ("c_matrix_dot_matrix_symm_upper", "matrix *a is empty.");
	if (b->size2 != a->size1)  c_error ("c_matrix_dot_matrix_symm_upper", "matrix size does not match.");

	return c_matrix_dsymm ("R", "U", alpha, b, a);
}

c_matrix *
c_matrix_dot_matrix_symm_lower (const double alpha, const c_matrix *b, const c_matrix *a)
{
	if (c_matrix_is_empty (b)) c_error ("c_matrix_dot_matrix_symm_lower", "matrix *b is empty.");
	if (c_matrix_is_empty (a)) c_error ("c_matrix_dot_matrix_symm_lower", "matrix *a is empty.");
	if (b->size2 != a->size1)  c_error ("c_matrix_dot_matrix_symm_lower", "matrix size does not match.");

	return c_matrix_dsymm ("R", "L", alpha, a, b);
}
