/*
 * c_matrix.c
 *
 *  Created on: 2014/04/03
 *      Author: utsugi
 */

#include <c_matrix.h>

/* blas */
extern void	dcopy_ (long *n, double *x, long *incx, double *y, long *incy);
extern void	c_error (const char * function_name, const char *error_msg);

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
c_matrix_alloc (const size_t size1, const size_t size2)
{
	c_matrix	*a;

	if (size1 <= 0 || size2 <= 0) c_error ("c_matrix_alloc", "invalid size.");

	a = _allocate_c_matrix ();
	if (!a) {
		fprintf (stderr, "WARNING: c_matrix_alloc: failed to allocate memory.\n");
		return NULL;
	}
	a->data = (double *) malloc (size1 * size2 * sizeof (double));
	if (!a->data) fprintf (stderr, "WARNING: c_matrix_alloc: failed to allocate memory.\n");
	else {
		a->size1 = size1;
		a->size2 = size2;
		a->lda = size1;
		a->tsize = a->lda * a->size2;
		a->owner = true;
	}
	return a;
}

c_matrix *
c_matrix_view_array (const size_t size1, const size_t size2, const size_t lda, double *data)
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
c_matrix_set (c_matrix *a, const int i, const int j, double val)
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
c_matrix_get_row (c_vector *y, const c_matrix *a, const size_t index)
{
	long	n;
	long	incx;
	long	incy;
	if (c_vector_is_empty (y)) c_error ("c_matrix_get_row", "vector is empty.");
	if (c_matrix_is_empty (a)) c_error ("c_matrix_get_row", "matrix is empty.");
	if (y->size != a->size1) c_error ("c_matrix_get_row", "vector and matrix size does not match.");
	if (index < 0 || index >= a->size1) c_error ("c_matrix_get_row", "index must be in [0, a->size1).");
	n = (long) a->size1;
	incx = (long) a->lda;
	incy = (long) y->stride;
	dcopy_ (&n, a->data + index, &incx, y->data, &incy);
	return;
}

/* copy column */
void
c_matrix_get_col (c_vector *y, const c_matrix *a, const size_t index)
{
	long	n;
	long	incx;
	long	incy;
	if (c_vector_is_empty (y)) c_error ("c_matrix_get_col", "vector is empty.");
	if (c_matrix_is_empty (a)) c_error ("c_matrix_get_col", "matrix is empty.");
	if (y->size != a->size1) c_error ("c_matrix_get_col", "vector and matrix size does not match.");
	if (index < 0 || index >= a->size2) c_error ("c_matrix_get_col", "index must be in [0, a->size2).");
	n = (long) a->size1;
	incx = 1;
	incy = (long) y->stride;
	dcopy_ (&n, a->data + a->lda * index, &incx, y->data, &incy);
	return;
}

void
c_matrix_set_row (c_matrix *a, const size_t index, const c_vector *x)
{
	long	n;
	long	incx;
	long	incy;
	if (c_vector_is_empty (x)) c_error ("c_matrix_set_row", "vector is empty.");
	if (c_matrix_is_empty (a)) c_error ("c_matrix_set_row", "matrix is empty.");
	if (x->size != a->size1) c_error ("c_matrix_set_row", "vector and matrix size does not match.");
	if (index < 0 || index >= a->size1) c_error ("c_matrix_set_row", "index must be in [0, a->size1).");
	n = (long) x->size;
	incx = (long) x->stride;
	incy = (long) a->lda;
	dcopy_ (&n, x->data, &incx, a->data + index, &incy);
	return;
}

void
c_matrix_set_col (c_matrix *a, const size_t index, const c_vector *x)
{
	long	n;
	long	incx;
	long	incy;
	if (c_vector_is_empty (x)) c_error ("c_matrix_set_col", "vector is empty.");
	if (c_matrix_is_empty (a)) c_error ("c_matrix_set_col", "matrix is empty.");
	if (x->size != a->size1) c_error ("c_matrix_set_col", "vector and matrix size does not match.");
	if (index < 0 || index >= a->size2) c_error ("c_matrix_set_col", "index must be in [0, a->size2).");
	n = (long) x->size;
	incx = (long) x->stride;
	incy = 1;
	dcopy_ (&n, x->data, &incx, a->data + a->lda * index, &incy);
	return;
}

/* create vector view of row */
c_vector *
c_matrix_row (c_matrix *a, int index)
{
	c_vector	*x;
	if (c_matrix_is_empty (a)) c_error ("c_matrix_row", "matrix is empty.");
	if (index < 0 || a->size1 <= index) c_error ("c_matrix_row", "index is invalid.");

	x = c_vector_view_array (a->size2, a->lda, a->data + index);
	return x;
}

/* create vector view of column */
c_vector *
c_matrix_column (c_matrix *a, int index)
{
	c_vector	*x;
	if (c_matrix_is_empty (a)) c_error ("c_matrix_column", "matrix is empty.");
	if (index < 0 || a->size2 <= index) c_error ("c_matrix_column", "index is invalid.");

	x = c_vector_view_array (a->size1, 1, a->data + index * a->lda);
	return x;
}

void
c_matrix_add_row (c_matrix *a)
{
	long	n;
	long	incx = 1;
	long	incy = 1;
	size_t	lda;

	if (c_matrix_is_empty (a)) c_error ("c_matrix_add_row", "matrix is empty.");

	lda = a->lda;
	a->size1++;
	a->lda++;
	a->tsize = a->lda * a->size2;
	a->data = (double *) realloc (a->data, a->tsize * sizeof (double));

	{
		int			j;
		c_vector	*col = c_vector_alloc (a->size1);

		n = (long) a->size1;
		for (j = a->size2 - 1; 0 < j; j--) {
			dcopy_ (&n, a->data + j * lda, &incx, col->data, &incy);
			dcopy_ (&n, col->data, &incy, a->data + j * a->lda, &incx);
		}
		c_vector_free (col);
		for (j = 0; j < a->size2; j++) c_matrix_set (a, a->size1 - 1, j, 0.);
	}
	return;
}

void
c_matrix_add_col (c_matrix *a)
{
	size_t	lda;

	if (c_matrix_is_empty (a)) c_error ("c_matrix_add_col", "matrix is empty.");

	lda = a->lda;
	a->size2++;
	a->tsize = a->lda * a->size2;
	a->data = (double *) realloc (a->data, a->tsize * sizeof (double));

	return;
}

void
c_matrix_remove_row (c_matrix *a)
{
	long	n;
	long	incx = 1;
	long	incy = 1;
	size_t	lda;

	if (c_matrix_is_empty (a)) c_error ("c_matrix_remove_row", "matrix is empty.");

	lda = a->lda;
	a->size1--;
	a->lda--;
	a->tsize = a->lda * a->size2;
	n = (long) a->size2;
	if (a->size1 > 0) {
		int			j;
		c_vector	*col = c_vector_alloc (a->size1);
		for (j = 1; j < a->size2; j++) {
			dcopy_ (&n, a->data + j * lda, &incx, col->data, &incy);
			dcopy_ (&n, col->data, &incx, a->data + j * a->lda, &incx);
		}
		c_vector_free (col);
	}
	if (a->data) a->data = (double *) realloc (a->data, a->tsize * sizeof (double));

	return;
}

void
c_matrix_remove_col (c_matrix *a)
{
	size_t	lda;

	if (c_matrix_is_empty (a)) c_error ("c_matrix_remove_col", "matrix is empty.");

	lda = a->lda;
	a->size2--;
	a->tsize = a->lda * a->size2;
	if (a->data) a->data = (double *) realloc (a->data, a->tsize * sizeof (double));

	return;
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
c_matrix_set_zero (c_matrix *a)
{
	int		i;
	if (c_matrix_is_empty (a)) c_error ("c_matrix_set_zero", "matrix is empty.");
	for (i = 0; i < a->tsize; i++) a->data[i] = 0.0;
	return;
}


c_matrix *
c_matrix_submatrix (const size_t size1, const size_t size2, const c_matrix *a)
{
	if (size1 < 0 || a->size1 < size1) c_error ("c_matrix_submatrix", "size1 must be in [0, a->size1]");
	if (size2 < 0 || a->size2 < size2) c_error ("c_matrix_submatrix", "size2 must be in [0, a->size2]");
	return c_matrix_view_array (size1, size2, a->lda, a->data);
}

void
c_matrix_fprintf (FILE *stream, const c_matrix *a, const char *format)
{
	int	i, j, k;
	if (c_matrix_is_empty (a)) c_error ("c_matrix_fprintf", "matrix is empty.");
	k = 0;
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
	int	i, j;
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
c_matrix_fread (FILE *stream, const size_t size1, const size_t size2)
{
	int			i;
	c_matrix	*a;
	a = c_matrix_alloc (size1, size2);
	for (i = 0; i < a->tsize; i++) {
		fscanf (stream, "%lf", &a->data[i]);
	}
	return a;
}
