/*
 * test_utils.c
 *
 *  Created on: 2014/04/14
 *      Author: utsugi
 */

#include <c_linalg.h>

/* lapack */
extern void	dlaruv_ (long seed[], long *n, double *x);

c_vector *
random_vector (size_t size)
{
	long		seed[4];
	long		n = (long) size;
	c_vector	*v = c_vector_alloc (size);
	dlaruv_ (seed, &n, v->data);
	return v;
}

c_matrix *
random_matrix (size_t size1, size_t size2)
{
	long		seed[4];
	long		n = (long) size1 * size2;
	c_matrix	*a = c_matrix_alloc (size1, size2);
	dlaruv_ (seed, &n, a->data);
	return a;
}

