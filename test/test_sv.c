/*
 * test_sv.c
 *
 *  Created on: 2014/04/14
 *      Author: utsugi
 */

#include <clinalg.h>

#include "test_clinalg.h"

extern size_t		size1;
extern size_t		size2;

/* check |a - u * diag(s) * vt| < 1.e-8 */
bool
test_SV_decomp (void)
{
	int			i;
	c_matrix	*a;
	c_vector	*s;
	c_matrix	*u;
	c_matrix	*vt;
	c_matrix	*c;

	double		nrm;

	a = random_matrix (size1, size2);

	/* svd */
	{
		c_matrix	*tmp = c_matrix_alloc (a->size1, a->size2);
		c_matrix_memcpy (tmp, a);
		c_linalg_SV_decomp (tmp, &u, &vt, &s, false);
		c_matrix_free (tmp);
	}

	/* u * diag(s) * vt */
	{
		c_matrix	*c1;
		c_matrix	*d = c_matrix_alloc (u->size2, vt->size1);
		c_vector	*diag = c_matrix_diagonal_view_array (d);
		c_matrix_set_zero (d);
		c_vector_memcpy (diag, s);
		c_vector_free (diag);
		c_vector_free (s);

		c1 = c_matrix_dot_matrix (1., u, d);
		c_matrix_free (u);
		c_matrix_free (d);

		c = c_matrix_dot_matrix (1., c1, vt);
		c_matrix_free (c1);
		c_matrix_free (vt);
	}

	c_matrix_sub (a, c);
	c_matrix_free (c);

	nrm = c_matrix_nrm (a, '1');
	c_matrix_free (a);

	return (nrm < 1.e-8);
}

bool
test_SV_solve (void)
{
	c_matrix	*a;
	c_vector	*x;
	c_vector	*y;
	c_vector	*b;
	double		nrm;

	a = random_matrix (size1, size2);
	x = random_vector (size2);
	y = c_matrix_dot_vector (1., a, x);
	c_vector_free (x);
	{
		int			rank;
		c_vector	*s;
		c_matrix	*tmp = c_matrix_alloc (a->size1, a->size2);
		c_matrix_memcpy (tmp, a);
		x = c_vector_alloc (y->size);
		c_vector_memcpy (x, y);
		c_linalg_SV_solve (1.e-8, tmp, x, &s, &rank);
		c_matrix_free (tmp);
	}
	b = c_matrix_dot_vector (1., a, x);
	c_matrix_free (a);
	c_vector_free (x);

	c_vector_axpy (-1., b, y);
	c_vector_free (b);

	nrm = c_vector_nrm (y);
	c_vector_free (y);

	return (nrm < 1.e-8);
}

bool
test_SV_lsd_solve (void)
{
	c_matrix	*a;
	c_vector	*x;
	c_vector	*y;
	c_vector	*b;
	double		nrm;

	a = random_matrix (size1, size2);
	x = random_vector (size2);
	y = c_matrix_dot_vector (1., a, x);
	c_vector_free (x);
	{
		int			rank;
		c_vector	*s;
		c_matrix	*tmp = c_matrix_alloc (a->size1, a->size2);
		c_matrix_memcpy (tmp, a);
		x = c_vector_alloc (y->size);
		c_vector_memcpy (x, y);
		c_linalg_SV_lsd_solve (1.e-8, tmp, x, &s, &rank);
		c_matrix_free (tmp);
	}
	b = c_matrix_dot_vector (1., a, x);
	c_matrix_free (a);
	c_vector_free (x);

	c_vector_axpy (-1., b, y);
	c_vector_free (b);

	nrm = c_vector_nrm (y);
	c_vector_free (y);

	return (nrm < 1.e-8);
}
