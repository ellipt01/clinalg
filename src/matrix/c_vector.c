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
	c_vector	*v = (c_vector *) malloc (sizeof (c_vector));
	if (!v) return NULL;
	v->size = 0;
	v->stride = 0;
	v->tsize = 0;
	v->data = NULL;
	v->owner = false;
	return v;
}

c_vector *
c_vector_alloc (const size_t size)
{
	c_vector	*v;

	if (size <= 0) c_error ("c_vector_alloc", "invalid size.");

	v = _allocate_c_vector ();
	if (!v) {
		fprintf (stderr, "WARNING: c_vector_alloc: failed to allocate memory.\n");
		return NULL;
	}
	v->data = (double *) malloc (size * sizeof (double));
	if (!v->data) fprintf (stderr, "WARNING: c_vector_alloc: failed to allocate memory.\n");
	else {
		v->size = size;
		v->stride = 1;
		v->tsize = v->stride * v->size;
		v->owner = true;
	}
	return v;
}

c_vector *
c_vector_view_array (const size_t size, const size_t stride, double *data)
{
	c_vector *v;

	if (!data) c_error ("c_vector_view_array", "array is empty.");
	if (size <= 0 || stride <= 0) c_error ("c_vector_view_array", "invalid size.");

	v = _allocate_c_vector ();
	v->size = size;
	v->stride = stride;
	v->tsize = v->stride * v->size;
	v->data = data;
	v->owner = false;
	return v;
}

bool
c_vector_is_empty (const c_vector *v)
{
	if (!v) return true;
	if (!v->data) return true;
	return false;
}

void
c_vector_free (c_vector *v)
{
	if (v) {
		if (v->data && v->owner) free (v->data);
		free (v);
	}
	return;
}

void
c_vector_set (c_vector *v, const int i, double val)
{
	int		index;
	if (c_vector_is_empty (v)) c_error ("c_vector_set", "vector is empty.");
	index = GET_INDEX_OF_VECTOR (v, i);
	if (index < 0 || v->size * v->stride <= index) c_error ("c_vector_set", "index out of range.");
	v->data[index] = val;
	return;
}

double
c_vector_get (const c_vector *v, const int i)
{
	int		index;
	if (c_vector_is_empty (v)) c_error ("c_vector_get", "vector is empty.");
	index = GET_INDEX_OF_VECTOR (v, i);
	if (index < 0 || v->tsize <= index) c_error ("c_vector_get", "index out of range.");
	return v->data[index];
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
c_vector_set_zero (c_vector *v)
{
	int		i;
	if (c_vector_is_empty (v)) c_error ("c_vector_set_zero", "vector is empty.");
	for (i = 0; i < v->size; i++) v->data[GET_INDEX_OF_VECTOR (v, i)] = 0.0;
	return;
}

void
c_vector_fprintf (FILE *stream, const c_vector *v, const char *format)
{
	int i;
	if (c_vector_is_empty (v)) c_error ("c_vector_fprintf", "vector is empty.");
	for (i = 0; i < v->size; i++) {
		fprintf (stream, format, v->data[GET_INDEX_OF_VECTOR (v, i)]);
		fprintf (stream, "\n");
	}
	return;
}
