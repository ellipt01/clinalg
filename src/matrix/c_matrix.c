/*
 * c_matrix.c
 *
 *  Created on: 2014/04/03
 *      Author: utsugi
 */

#include <c_matrix.h>
#include "../../include/clinalg_macros.h"
#include "../../include/clinalg_utils.h"

#include "private.h"

static c_matrix *
_allocate_c_matrix (void)
{
	c_matrix	*a = (c_matrix *) malloc (sizeof (c_matrix));
	if (!a) return NULL;
	a->size1 = 0;
	a->size2 = 0;
	a->lda = 0;
	a->tsize = 0;
	a->data = NULL;
	a->owner = false;
	return a;
}

c_matrix *
c_matrix_alloc (const int size1, const int size2)
{
	c_matrix	*a;

	if (size1 <= 0 || size2 <= 0) c_error ("c_matrix_alloc", "invalid size.");

	a = _allocate_c_matrix ();
	if (!a) c_error ("c_matrix_alloc", "failed to allocate memory.");
	a->data = (double *) malloc (size1 * size2 * sizeof (double));
	if (!a->data) c_error ("c_matrix_alloc", "failed to allocate memory.");
	a->size1 = size1;
	a->size2 = size2;
	a->lda = size1;
	a->tsize = a->lda * a->size2;
	a->owner = true;

	return a;
}

void
c_matrix_realloc (const int tsize, c_matrix *x, const int size1, const int size2)
{
	if (c_matrix_is_empty (x)) c_error ("c_matrix_realloc", "matrix is empty.");
	if (!x->owner) c_error ("c_matrix_realloc", "cannot reallocate matrix of !x->owner.");
	if (x->tsize == tsize) return;

	x->data = (double *) realloc (x->data, tsize * sizeof (double));
	if (x->data == NULL) c_error ("c_matrix_realloc", "reallocation of array failed.");
	x->tsize = tsize;
	if (x->size1 != size1) x->size1 = size1;
	if (x->size2 != size2) x->size2 = size2;

	return;
}

c_matrix *
c_matrix_view_array (const int size1, const int size2, const int lda, double *data)
{
	c_matrix *a;

	if (!data) c_error ("c_matrix_view_array", "array is empty.");
	if (size1 <= 0 || size2 <= 0 || lda <= 0) c_error ("c_matrix_view_array", "invalid size.");

	a = _allocate_c_matrix ();
	a->size1 = size1;
	a->size2 = size2;
	a->lda = lda;
	a->tsize = a->lda * a->size2;
	a->data = data;
	a->owner = false;
	return a;
}

bool
c_matrix_is_empty (const c_matrix *a)
{
	if (!a) return true;
	if (!a->data) return true;
	if (a->size1 <= 0 || a->size2 <= 0) return true;
	return false;
}

bool
c_matrix_is_square (const c_matrix *a)
{
	return (a->size1 == a->size2);
}

void
c_matrix_free (c_matrix *a)
{
	if (a) {
		if (a->data && a->owner) free (a->data);
		free (a);
	}
	return;
}

void
c_matrix_set (c_matrix *a, const int i, const int j, const double val)
{
	int		index;
	if (c_matrix_is_empty (a)) c_error ("c_matrix_set", "matrix is empty.");
	index = INDEX_OF_MATRIX (a, i, j);
	if (index < 0 || a->tsize <= index) c_error ("c_matrix_set", "index out of range.");
	a->data[index] = val;
	return;
}

double
c_matrix_get (const c_matrix *a, const int i, const int j)
{
	int		index;
	if (c_matrix_is_empty (a)) c_error ("c_matrix_get", "matrix is empty.");
	index = INDEX_OF_MATRIX (a, i, j);
	if (index < 0 || a->tsize <= index) c_error ("c_matrix_get", "index out of range.");
	return a->data[index];
}

/* copy row */
void
c_matrix_get_row (c_vector *y, const c_matrix *a, const int index)
{
	int		n;
	int		incx;
	int		incy;
	if (c_vector_is_empty (y)) c_error ("c_matrix_get_row", "vector is empty.");
	if (c_matrix_is_empty (a)) c_error ("c_matrix_get_row", "matrix is empty.");
	if (y->size != a->size2) c_error ("c_matrix_get_row", "vector and matrix size does not match.");
	if (index < 0 || index >= a->size1) c_error ("c_matrix_get_row", "index must be in [0, a->size1).");
	n = (int) a->size2;
	incx = (int) a->lda;
	incy = (int) y->stride;
	F77CALL (dcopy) (&n, a->data + index, &incx, y->data, &incy);
	return;
}

/* copy column */
void
c_matrix_get_col (c_vector *y, const c_matrix *a, const int index)
{
	int		n;
	int		incx;
	int		incy;
	if (c_vector_is_empty (y)) c_error ("c_matrix_get_col", "vector is empty.");
	if (c_matrix_is_empty (a)) c_error ("c_matrix_get_col", "matrix is empty.");
	if (y->size != a->size1) c_error ("c_matrix_get_col", "vector and matrix size does not match.");
	if (index < 0 || index >= a->size2) c_error ("c_matrix_get_col", "index must be in [0, a->size2).");
	n = (int) a->size1;
	incx = 1;
	incy = (int) y->stride;
	F77CALL (dcopy) (&n, a->data + a->lda * index, &incx, y->data, &incy);
	return;
}

void
c_matrix_set_row (c_matrix *a, const int index, const c_vector *x)
{
	int		n;
	int		incx;
	int		incy;
	if (c_vector_is_empty (x)) c_error ("c_matrix_set_row", "vector is empty.");
	if (c_matrix_is_empty (a)) c_error ("c_matrix_set_row", "matrix is empty.");
	if (x->size != a->size2) c_error ("c_matrix_set_row", "vector and matrix size does not match.");
	if (index < 0 || index >= a->size1) c_error ("c_matrix_set_row", "index must be in [0, a->size1).");
	n = (int) x->size;
	incx = (int) x->stride;
	incy = (int) a->lda;
	F77CALL (dcopy) (&n, x->data, &incx, a->data + index, &incy);
	return;
}

void
c_matrix_set_col (c_matrix *a, const int index, const c_vector *x)
{
	int		n;
	int		incx;
	int		incy;
	if (c_vector_is_empty (x)) c_error ("c_matrix_set_col", "vector is empty.");
	if (c_matrix_is_empty (a)) c_error ("c_matrix_set_col", "matrix is empty.");
	if (x->size != a->size1) c_error ("c_matrix_set_col", "vector and matrix size does not match.");
	if (index < 0 || index >= a->size2) c_error ("c_matrix_set_col", "index must be in [0, a->size2).");
	n = (int) x->size;
	incx = (int) x->stride;
	incy = 1;
	F77CALL (dcopy) (&n, x->data, &incx, a->data + a->lda * index, &incy);
	return;
}

/* create vector view of row */
c_vector *
c_matrix_row (c_matrix *a, const int index)
{
	c_vector	*x;
	if (c_matrix_is_empty (a)) c_error ("c_matrix_row", "matrix is empty.");
	if (index < 0 || a->size1 <= index) c_error ("c_matrix_row", "index is invalid.");

	x = c_vector_view_array (a->size2, a->lda, a->data + index);
	return x;
}

/* create vector view of column */
c_vector *
c_matrix_column (c_matrix *a, const int index)
{
	c_vector	*x;
	if (c_matrix_is_empty (a)) c_error ("c_matrix_column", "matrix is empty.");
	if (index < 0 || a->size2 <= index) c_error ("c_matrix_column", "index is invalid.");

	x = c_vector_view_array (a->size1, 1, a->data + index * a->lda);
	return x;
}


void
c_matrix_memcpy (c_matrix *dest, const c_matrix *src)
{
	int			j;
	c_vector	*col;
	if (c_matrix_is_empty (src)) c_error ("c_matrix_memcpy", "first matrix is empty.");
	if (c_matrix_is_empty (dest)) c_error ("c_matrix_memcpy", "second matrix is empty.");
	if (dest->size1 != src->size1 || dest->size2 != src->size2) c_error ("c_matrix_memcpy", "matrix size does not match.");
	col = c_vector_alloc (src->size1);
	for (j = 0; j < src->size2; j++) {
		c_matrix_get_col (col, src, j);
		c_matrix_set_col (dest, j, col);
	}
	c_vector_free (col);
	return;
}

void
c_matrix_mncopy (c_matrix *dest, const int m0, const int n0, const int m, const int n, const c_matrix *src)
{
	int		j;
	int		len;
	int		inc = 1;
	if (m0 < 0 || src->size1 - 1 <= m0) c_error ("c_matrix_mncopy", "m0 must be in [0, src->size1 - 1)");
	if (n0 < 0 || src->size2 - 1 <= n0) c_error ("c_matrix_mncopy", "n0 must be in [0, src->size2 - 1)");
	if (m0 + m < 0 || src->size1 < m0 + m) c_error ("c_matrix_mncopy", "m0 + m must be in [0, src->size1]");
	if (n0 + n < 0 || src->size2 < n0 + n) c_error ("c_matrix_mncopy", "n0 + n must be in [0, src->size2]");
	if (dest->size1 < m || dest->size2 < n) c_error ("c_matrix_mncopy", "index out of range.");

	len = (int) m;
	for (j = 0; j < n; j++) {
		F77CALL (dcopy) (&len, POINTER_OF_MATRIX (src, m0, n0 + j), &inc, POINTER_OF_MATRIX (dest, 0, j), &inc);
	}
	return;
}

void
c_matrix_upper_triangular_memcpy (c_matrix *tr, const c_matrix *a)
{
	int		j;
	int		incx = 1;
	int		incy = 1;

	int	min_m = C_MIN (tr->size1, a->size1);
	int	min_n = C_MIN (tr->size2, a->size2);
	for (j = 0; j < min_n; j++) {
		int	n = (j + 1 < min_m) ? (int) (j + 1) : (int) min_m;
		F77CALL (dcopy) (&n, a->data + j * a->lda, &incx, tr->data + j * tr->lda, &incy);
	}
	return;
}

void
c_matrix_lower_triangular_memcpy (c_matrix *tr, const c_matrix *a)
{
	int		j;
	int		incx = 1;
	int		incy = 1;

	int	min_m = C_MIN (tr->size1, a->size1);
	int	min_n = C_MIN (tr->size2, a->size2);
	for (j = 0; j < min_n; j++) {
		int	n;
		n = min_m - j;
		if (n <= 0) break;
		F77CALL (dcopy) (&n, a->data + j * (a->lda + 1), &incx, tr->data + j * (tr->lda + 1), &incy);
	}
	return;
}

void
c_matrix_set_all (c_matrix *a, const double val)
{
	int			i, j;
	if (c_matrix_is_empty (a)) c_error ("c_matrix_set_all", "matrix is empty.");
	for (i = 0; i < a->size1; i++) {
		for (j = 0; j < a->size2; j++) c_matrix_set (a, i, j, val);
	}
	return;
}

void
c_matrix_set_zero (c_matrix *a)
{
	if (c_matrix_is_empty (a)) c_error ("c_matrix_set_zero", "matrix is empty.");
	c_matrix_set_all (a, 0.);
	return;
}

c_vector *
c_matrix_get_diagonal (const c_matrix *a)
{
	int			n;
	int			inc;
	int			lda;
	c_vector	*d;

	if (c_matrix_is_empty (a)) c_error ("c_matrix_get_diagonal", "matrix is empty.");

	n = (int) C_MIN (a->size1, a->size2);
	lda = (int) a->lda + 1;
	d = c_vector_alloc (n);
	inc = 1;
	F77CALL (dcopy) (&n, a->data, &lda, d->data, &inc);

	return d;
}

c_vector *
c_matrix_diagonal_view_array (const c_matrix *a)
{
	int		min_mn;
	c_vector	*d;

	if (c_matrix_is_empty (a)) c_error ("c_matrix_get_diagonal_view_array", "matrix is empty.");

	min_mn = C_MIN (a->size1, a->size2);
	d = c_vector_view_array (min_mn, a->lda + 1, a->data);

	return d;
}

void
c_matrix_set_diagonal (const c_vector *d, c_matrix *a)
{
	int			n;
	int			inc;
	int			lda;

	if (c_matrix_is_empty (a)) c_error ("c_matrix_set_diagonal", "matrix is empty.");

	n = (int) C_MIN (a->size1, a->size2);
	inc = (int) d->stride;
	lda = (int) a->lda + 1;
	F77CALL (dcopy) (&n, d->data, &inc, a->data, &lda);

	return;
}

c_matrix *
c_matrix_submatrix (const int m0, const int n0, const int m, const int n, const c_matrix *a)
{
	if (m0 < 0 || a->size1 - 1 < m0) c_error ("c_matrix_submatrix", "m0 must be in [0, a->size1 - 1)");
	if (n0 < 0 || a->size2 - 1 < n0) c_error ("c_matrix_submatrix", "n0 must be in [0, a->size2 - 1)");
	if (m0 + m < 0 || a->size1 < m0 + m) c_error ("c_matrix_submatrix", "m0 + m must be in [0, a->size1]");
	if (n0 + n < 0 || a->size2 < n0 + n) c_error ("c_matrix_submatrix", "n0 + n must be in [0, a->size2]");
	return c_matrix_view_array (m, n, a->lda, &a->data[INDEX_OF_MATRIX(a, m0, n0)]);
}

void
c_matrix_fprintf (FILE *stream, const c_matrix *a, const char *format)
{
	int		i, j;
	if (c_matrix_is_empty (a)) c_error ("c_matrix_fprintf", "matrix is empty.");
	for (j = 0; j < a->size2; j++) {
		for (i = 0; i < a->size1; i++) {
			fprintf (stream, format, c_matrix_get (a, i, j));
			fprintf (stream, "\n");
		}
	}
	return;
}

void
c_matrix_fprintf2 (FILE *stream, const c_matrix *a, const char *format)
{
	int		i, j;
	if (c_matrix_is_empty (a)) c_error ("c_matrix_fprintf2", "matrix is empty.");
	for (i = 0; i < a->size1; i++) {
		for (j = 0; j < a->size2; j++) {
			fprintf (stream, format, c_matrix_get (a, i, j));
			if (j < a->size2 - 1) fprintf (stream, " ");
		}
		fprintf (stream, "\n");
	}
	return;
}

c_matrix *
c_matrix_fread (FILE *stream, const int size1, const int size2)
{
	int			i, k;
	c_matrix	*a;
	a = c_matrix_alloc (size1, size2);
	k = 0;
	for (i = 0; i < a->tsize; i++) {
		if (fscanf (stream, "%lf", &a->data[k]) < 0) continue;
		k++;
	}
	return a;
}
