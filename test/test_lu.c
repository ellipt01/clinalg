/*
 * test_lu.c
 *
 *  Created on: 2014/04/17
 *      Author: utsugi
 */

#include <c_linalg.h>

#include "test_clinalg.h"

bool
test_LU_decomp (void)
{
	int				i;
	size_t			size = 50;
	c_matrix		*a;
	c_matrix		*lu;
	c_matrix		*l;
	c_matrix		*u;
	c_matrix		*b;
	c_matrix		*ap;
	c_vector_int	*p;
	c_matrix		*perm;
	double			nrm;

	a = random_matrix (size, size);

	lu = c_matrix_alloc (a->size1, a->size2);
	c_matrix_memcpy (lu, a);
	c_linalg_LU_decomp (lu, &p);

	l = c_matrix_alloc (size, size);
	c_matrix_set_zero (l);
	c_matrix_lower_triangular_memcpy (l, lu);
	for (i = 0; i < size; i++) c_matrix_set (l, i, i, 1.);

	u = c_matrix_alloc (size, size);
	c_matrix_set_zero (u);
	c_matrix_upper_triangular_memcpy (u, lu);
	c_matrix_free (lu);

	perm = c_matrix_identity (a->size1, a->size2);
	for (i = 0; i < size; i++) {
		int		pi = (size_t) c_vector_int_get (p, i) - 1;
		if (pi != i) c_matrix_swap_rows (i, pi, perm);
	}
	c_vector_int_free (p);

	b = c_matrix_dot_matrix (1., l, u, 0.);
	c_matrix_free (l);
	c_matrix_free (u);

	ap = c_matrix_dot_matrix (1., perm, a, 0.);
	c_matrix_free (a);
	c_matrix_free (perm);

	c_matrix_sub (b, ap);
	c_matrix_free (ap);

	nrm = c_matrix_nrm (b, '1');
	c_matrix_free (b);

	return (nrm < 1.e-8);
}

bool
test_LU_solve (void)
{
	size_t			size = 50;
	c_matrix		*a;
	c_vector		*x;
	c_vector		*b;

	c_matrix		*lu;
	c_vector		*y;
	c_vector_int	*p;
	double			nrm;

	a = random_matrix (size, size);
	x = random_vector (size);
	b = c_matrix_dot_vector (1., a, x, 0.);
	c_vector_free (x);

	lu = c_matrix_alloc (a->size1, a->size2);
	c_matrix_memcpy (lu, a);
	c_linalg_LU_decomp (lu, &p);

	x = c_vector_alloc (b->size);
	c_vector_memcpy (x, b);
	c_linalg_LU_solve (lu, x, p);
	c_matrix_free (lu);
	c_vector_int_free (p);

	y = c_matrix_dot_vector (1., a, x, 0.);
	c_matrix_free (a);
	c_vector_free (x);

	c_vector_sub (b, y);
	c_vector_free (y);

	nrm = c_vector_nrm (b);
	c_vector_free (b);

	return (nrm < 1.e-8);
}
