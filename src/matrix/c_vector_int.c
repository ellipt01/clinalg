/*
 * c_vector_int.c
 *
 *  Created on: 2014/04/03
 *      Author: utsugi
 */

#include <c_matrix.h>

extern void	c_error (const char * function_name, const char *error_msg);

static c_vector_int *
_allocate_c_vector_int (void)
{
	c_vector_int	*v = (c_vector_int *) malloc (sizeof (c_vector_int));
	if (!v) return NULL;
	v->size = 0;
	v->stride = 0;
	v->tsize = 0;
	v->data = NULL;
	v->owner = false;
	return v;
}

c_vector_int *
c_vector_int_alloc (const size_t size)
{
	c_vector_int	*v;

	if (size <= 0) c_error ("c_vector_int_alloc", "invalid size.");

	v = _allocate_c_vector_int ();
	if (!v) {
		fprintf (stderr, "WARNING: c_vector_int_alloc: failed to allocate memory.\n");
		return NULL;
	}
	v->data = (int *) malloc (size * sizeof (int));
	if (!v->data) fprintf (stderr, "WARNING: c_vector_int_alloc: failed to allocate memory.\n");
	else {
		v->size = size;
		v->stride = 1;
		v->tsize = v->stride * v->size;
		v->owner = true;
	}
	return v;
}

c_vector_int *
c_vector_int_view_array (const size_t size, const size_t stride, int *data)
{
	c_vector_int *v;

	if (!data) c_error ("c_vector_int_view_array", "array is empty.");
	if (size <= 0 || stride <= 0) c_error ("c_vector_int_view_array", "invalid size.");

	v = _allocate_c_vector_int ();
	v->size = size;
	v->stride = stride;
	v->tsize = v->stride * v->size;
	v->data = data;
	v->owner = false;
	return v;
}

bool
c_vector_int_is_empty (const c_vector_int *v)
{
	if (!v) return true;
	if (!v->data) return true;
	return false;
}

void
c_vector_int_free (c_vector_int *v)
{
	if (v) {
		if (v->data && v->owner) free (v->data);
		free (v);
	}
	return;
}

void
c_vector_int_set (c_vector_int *v, const int i, int val)
{
	int		index;
	if (c_vector_int_is_empty (v)) c_error ("c_vector_int_set", "vector is empty.");
	index = INDEX_OF_VECTOR (v, i);
	if (index < 0 || v->size * v->stride <= index) c_error ("c_vector_int_set", "index out of range.");
	v->data[index] = val;
	return;
}

int
c_vector_int_get (const c_vector_int *v, const int i)
{
	int		index;
	if (c_vector_int_is_empty (v)) c_error ("c_vector_int_get", "vector is empty.");
	index = INDEX_OF_VECTOR (v, i);
	if (index < 0 || v->tsize <= index) c_error ("c_vector_int_get", "index out of range.");
	return v->data[index];
}

void
c_vector_int_memcpy (c_vector_int *dest, const c_vector_int *src)
{
	int		i;
	long	n;
	long	incx;
	long	incy;
	if (c_vector_int_is_empty (src)) c_error ("c_vector_int_memcpy", "first vector is empty.");
	if (c_vector_int_is_empty (dest)) c_error ("c_vector_int_memcpy", "second vector is empty.");
	if (dest->size != src->size) c_error ("c_vector_int_memcpy", "vector size does not match.");
	n = (long) src->size;
	incx = (long) src->stride;
	incy = (long) dest->stride;
	for (i = 0; i < src->size; i++) dest[i * src->stride] = src[i * src->stride];
	return;
}

void
c_vector_int_set_zero (c_vector_int *v)
{
	int		i;
	if (c_vector_int_is_empty (v)) c_error ("c_vector_int_set_zero", "vector is empty.");
	for (i = 0; i < v->size; i++) v->data[INDEX_OF_VECTOR (v, i)] = 0;
	return;
}

void
c_vector_int_fprintf (FILE *stream, const c_vector_int *v, const char *format)
{
	int i;
	if (c_vector_int_is_empty (v)) c_error ("c_vector_int_fprintf", "vector is empty.");
	for (i = 0; i < v->size; i++) {
		fprintf (stream, format, v->data[INDEX_OF_VECTOR (v, i)]);
		fprintf (stream, "\n");
	}
	return;
}
