/*
 * c_vector.c
 *
 *  Created on: 2014/04/03
 *      Author: utsugi
 */

#include <c_matrix.h>

/* blas */
extern void	dcopy_ (long *n, double *x, long *incx, double *y, long *incy);
extern void	c_error (const char * function_name, const char *error_msg);

static c_vector *
_allocate_c_vector (void)
{
	c_vector	*x = (c_vector *) malloc (sizeof (c_vector));
	if (!x) return NULL;
	x->size = 0;
	x->stride = 0;
	x->tsize = 0;
	x->data = NULL;
	x->owner = false;
	return x;
}

c_vector *
c_vector_alloc (const size_t size)
{
	c_vector	*x;

	if (size <= 0) c_error ("c_vector_alloc", "invalid size.");

	x = _allocate_c_vector ();
	if (!x) {
		fprintf (stderr, "WARNING: c_vector_alloc: failed to allocate memory.\n");
		return NULL;
	}
	x->data = (double *) malloc (size * sizeof (double));
	if (!x->data) fprintf (stderr, "WARNING: c_vector_alloc: failed to allocate memory.\n");
	else {
		x->size = size;
		x->stride = 1;
		x->tsize = x->stride * x->size;
		x->owner = true;
	}
	return x;
}

void
c_vector_realloc (c_vector *x, const size_t size)
{
	void	*data;
	size_t	memsize;

	if (c_vector_is_empty (x)) c_error ("c_vector_realloc", "vector is empty.");
	if (!x->owner) c_error ("c_vector_realloc", "cannot reallocate vector your own.");
	memsize = size * sizeof (double);
	if ((data = (double *) realloc (x->data, memsize)) == NULL) c_error ("c_vector_realloc", "cannot reallocate vector.");
	x->data = data;
	x->size = size;

	return;
}

c_vector *
c_vector_view_array (const size_t size, const size_t stride, double *data)
{
	c_vector *x;

	if (!data) c_error ("c_vector_view_array", "array is empty.");
	if (size <= 0 || stride <= 0) c_error ("c_vector_view_array", "invalid size.");

	x = _allocate_c_vector ();
	x->size = size;
	x->stride = stride;
	x->tsize = x->stride * x->size;
	x->data = data;
	x->owner = false;
	return x;
}

bool
c_vector_is_empty (const c_vector *x)
{
	if (!x) return true;
	if (!x->data) return true;
	return false;
}

void
c_vector_free (c_vector *x)
{
	if (x) {
		if (x->data && x->owner) free (x->data);
		free (x);
	}
	return;
}

void
c_vector_set (c_vector *x, const int i, double val)
{
	int		index;
	if (c_vector_is_empty (x)) c_error ("c_vector_set", "vector is empty.");
	index = INDEX_OF_VECTOR (x, i);
	if (index < 0 || x->size * x->stride <= index) c_error ("c_vector_set", "index out of range.");
	x->data[index] = val;
	return;
}

double
c_vector_get (const c_vector *x, const int i)
{
	int		index;
	if (c_vector_is_empty (x)) c_error ("c_vector_get", "vector is empty.");
	index = INDEX_OF_VECTOR (x, i);
	if (index < 0 || x->tsize <= index) c_error ("c_vector_get", "index out of range.");
	return x->data[index];
}

void
c_vector_memcpy (c_vector *dest, const c_vector *src)
{
	long	n;
	long	incx;
	long	incy;
	if (c_vector_is_empty (src)) c_error ("c_vector_memcpy", "first vector is empty.");
	if (c_vector_is_empty (dest)) c_error ("c_vector_memcpy", "second vector is empty.");
	if (dest->size != src->size) c_error ("c_vector_memcpy", "vector size does not match.");
	n = (long) src->size;
	incx = (long) src->stride;
	incy = (long) dest->stride;
	dcopy_ (&n, src->data, &incx, dest->data, &incy);
	return;
}

void
c_vector_set_zero (c_vector *x)
{
	int		i;
	if (c_vector_is_empty (x)) c_error ("c_vector_set_zero", "vector is empty.");
	for (i = 0; i < x->size; i++) x->data[INDEX_OF_VECTOR (x, i)] = 0.0;
	return;
}

c_vector *
c_vector_subvector (const size_t size, const c_vector *x)
{
	if (size < 0 || x->size < size) c_error ("c_vector_subvector", "size must be in [0, x->size]");
	return c_vector_view_array (size, x->stride, x->data);
}

void
c_vector_fprintf (FILE *stream, const c_vector *x, const char *format)
{
	int i;
	if (c_vector_is_empty (x)) c_error ("c_vector_fprintf", "vector is empty.");
	for (i = 0; i < x->size; i++) {
		fprintf (stream, format, x->data[INDEX_OF_VECTOR (x, i)]);
		fprintf (stream, "\n");
	}
	return;
}
