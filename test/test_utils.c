/*
 * test_utils.c
 *
 *  Created on: 2014/04/14
 *      Author: utsugi
 */

#include <clinalg.h>

c_vector *
random_vector (int size)
{
	int			i;
	c_vector	*v = c_vector_alloc (size);
	for (i = 0; i < size; i++) {
		double	r = (double) rand () / RAND_MAX;
		v->data[i] = 1. - 2. * r;
	}
	return v;
}

c_matrix *
random_matrix (int size1, int size2)
{
	int			i;
	c_matrix	*a = c_matrix_alloc (size1, size2);
	for (i = 0; i < size1 * size2; i++) {
		double	r = (double) rand () / RAND_MAX;
		a->data[i] = 1. - 2. * r;
	}
	return a;
}

