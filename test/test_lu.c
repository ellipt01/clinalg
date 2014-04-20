/*
 * test_lu.c
 *
 *  Created on: 2014/04/17
 *      Author: utsugi
 */

#include <math.h>
#include <c_linalg.h>

#include "test_clinalg.h"

bool
test_LU_decomp (void)
{
	int				i;
	size_t			size1 = 50;
	size_t			size2 = 60;
	c_matrix		*a;
	c_matrix		*lu;
	c_matrix		*l;
	c_matrix		*u;
	c_matrix		*b;
	c_matrix		*ap;
	c_vector_int	*p;
	c_matrix		*perm;
	double			nrm;

	a = random_matrix (size1, size2);

	lu = c_matrix_alloc (a->size1, a->size2);
	c_matrix_memcpy (lu, a);
	c_linalg_LU_decomp (lu, &p);

	c_linalg_LU_unpack (lu, &l, &u);
	c_matrix_free (lu);

	perm = c_linalg_permutation_matrix_row (a->size1, p);
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

bool
test_LU_invert (void)
{
	size_t			size = 50;
	c_matrix		*a;
	c_matrix		*lu;
	c_matrix		*c;
	c_vector_int	*p;
	double			nrm;

	a = random_matrix (size, size);

	lu = c_matrix_alloc (a->size1, a->size2);
	c_matrix_memcpy (lu, a);
	c_linalg_LU_decomp (lu, &p);
	c_linalg_LU_invert (lu, p);
	c_vector_int_free (p);

	c = c_matrix_dot_matrix (1., lu, a, 0.);
	c_matrix_free (a);
	c_matrix_free (lu);

	nrm = c_matrix_nrm (c, '1');
	c_matrix_free (c);

	return (fabs (nrm - 1.) < 1.e-8);
}

bool
test_LU_1up (void)
{
	double			nrm;
	size_t			size1 = 60;
	size_t			size2 = 50;
	c_matrix		*a;
	c_vector		*t;
	c_vector		*s;

	c_matrix		*l;
	c_matrix		*u;
	c_vector_int	*p;

	a = random_matrix (size1, size2);
	{
		c_matrix		*lu = c_matrix_alloc (a->size1, a->size2);
		c_matrix_memcpy (lu, a);
		c_linalg_LU_decomp (lu, &p);
		c_linalg_LU_unpack (lu, &l, &u);
		c_matrix_free (lu);
	}

	s = random_vector (size1);
	t = random_vector (size2);
	{
		c_matrix	*s0 = c_matrix_view_array (s->size, 1, s->size, s->data);
		c_matrix	*t0 = c_matrix_view_array (t->size, 1, t->size, t->data);
		c_matrix	*c = c_matrix_dot_matrix_transpose (1., s0, t0, 0.);
		c_matrix_free (s0);
		c_matrix_free (t0);
		c_matrix_add (a, c);
		c_matrix_free (c);
	}

	c_linalg_LU_1up (l, u, p, s, t);
	c_vector_free (s);
	c_vector_free (t);

	{
		c_matrix	*a1 = c_matrix_dot_matrix (1., l, u, 0.);
		c_matrix	*perm = c_linalg_permutation_matrix_row (a1->size1, p);
		c_matrix	*a2 = c_matrix_dot_matrix (1., perm, a1, 0.);
		c_matrix_free (a1);
		c_matrix_free (l);
		c_matrix_free (u);

		c_matrix_sub (a, a2);
		c_matrix_free (a2);
	}
	nrm = c_matrix_nrm (a, '1');
	c_matrix_free (a);

	return nrm;
}
