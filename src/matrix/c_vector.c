/*
 * c_vector.c
 *
 *  Created on: 2014/04/03
 *      Author: utsugi
 */

#include <c_linalg_macros.h>
#include <c_vector.h>
#include <c_linalg_utils.h>

#include "private.h"

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
c_vector_alloc (const int size)
{
	c_vector	*x;

	if (size <= 0) c_error ("c_vector_alloc", "invalid size.");

	x = _allocate_c_vector ();
	if (!x) c_error ("c_vector_alloc", "failed to allocate memory.");
	x->data = (double *) malloc (size * sizeof (double));
	if (!x->data) c_error ("c_vector_alloc", "failed to allocate memory.");
	x->size = size;
	x->stride = 1;
	x->tsize = x->stride * x->size;
	x->owner = true;

	return x;
}

void
c_vector_realloc (const int tsize, c_vector *x, const int size)
{
	if (c_vector_is_empty (x)) c_error ("c_vector_realloc", "vector is empty.");
	if (!x->owner) c_error ("c_vector_realloc", "cannot reallocate vector of !x->owner.");
	if (x->tsize == tsize) return;

	x->data = (double *) realloc (x->data, tsize * sizeof (double));
	if (x->data == NULL) c_error ("c_vector_realloc", "reallocation of array failed.");
	x->tsize = tsize;
	if (x->size != size) x->size = size;

	return;
}

c_vector *
c_vector_view_array (const int size, const int stride, double *data)
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
	if (x->size <= 0) return true;
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
c_vector_set (c_vector *x, const int i, const double val)
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
	if (c_vector_is_empty (src)) c_error ("c_vector_memcpy", "first vector is empty.");
	if (c_vector_is_empty (dest)) c_error ("c_vector_memcpy", "second vector is empty.");
	if (dest->size != src->size) c_error ("c_vector_memcpy", "vector size does not match.");
	dcopy_ (&src->size, src->data, &src->stride, dest->data, &dest->stride);
	return;
}

void
c_vector_ncopy (c_vector *dest, const int n0, const int n, const c_vector *src)
{
	if (c_vector_is_empty (src)) c_error ("c_vector_ncopy", "first vector is empty.");
	if (c_vector_is_empty (dest)) c_error ("c_vector_ncopy", "second vector is empty.");
	if (n0 < 0 || src->size - 1 <= n0) c_error ("c_vector_ncopy", "n0 must be in [0, size - 1).");
	if (n0 < 0 || src->size - 1 <= n0) c_error ("c_vector_ncopy", "n0 must be in [0, size - 1).");
	if (n0 + n < 0 || src->size < n0 + n) c_error ("c_vector_ncopy", "n0 + n must be in [0, src->size]");
	if (dest->size < n) c_error ("c_vector_ncopy", "index out of range.");

	dcopy_ (&n, src->data + n0, &src->stride, dest->data, &dest->stride);
	return;
}

void
c_vector_set_all (c_vector *x, const double val)
{
	int		i;
	if (c_vector_is_empty (x)) c_error ("c_vector_set_all", "first vector is empty.");
	for (i = 0; i < x->size; i++) x->data[INDEX_OF_VECTOR (x, i)] = val;
}

void
c_vector_set_zero (c_vector *x)
{
	if (c_vector_is_empty (x)) c_error ("c_vector_set_zero", "vector is empty.");
	c_vector_set_all (x, 0.);
	return;
}

c_vector *
c_vector_subvector (const int size, const c_vector *x)
{
	if (size < 0 || x->size < size) c_error ("c_vector_subvector", "size must be in [0, x->size]");
	return c_vector_view_array (size, x->stride, x->data);
}

void
c_vector_fprintf (FILE *stream, const c_vector *x, const char *format)
{
	int 	i;
	if (c_vector_is_empty (x)) c_error ("c_vector_fprintf", "vector is empty.");
	for (i = 0; i < x->size; i++) {
		fprintf (stream, format, x->data[INDEX_OF_VECTOR (x, i)]);
		fprintf (stream, "\n");
	}
	return;
}
