/*
 * test_lu.c
 *
 *  Created on: 2014/04/17
 *      Author: utsugi
 */

#include <math.h>
#include <clinalg.h>

#include "test_clinalg.h"

extern int		size1;
extern int		size2;

bool
test_LU_decomp (void)
{
	c_matrix		*a;
	c_matrix		*lu;
	c_matrix		*l;
	c_matrix		*u;
	c_matrix		*b;
	c_vector_int	*p;
	double			nrm;

	a = random_matrix (size1, size2);

	lu = c_matrix_alloc (a->size1, a->size2);
	c_matrix_memcpy (lu, a);
	c_linalg_LU_decomp (lu, &p);
	c_linalg_LU_unpack (lu, &l, &u);
	c_matrix_free (lu);

	b = c_matrix_dot_matrix (1., l, u);
	c_matrix_free (l);
	c_matrix_free (u);

	c_matrix_permute_rows (a, p);
	c_vector_int_free (p);

	c_matrix_sub (a, b);
	c_matrix_free (b);

	nrm = c_matrix_nrm (a, '1');
	c_matrix_free (a);

	return (nrm < 1.e-8);
}

bool
test_LU_solve (void)
{
	c_matrix		*a;
	c_vector		*x;
	c_vector		*b;

	c_matrix		*lu;
	c_vector		*y;
	double			nrm;

	a = random_matrix (size1, size1);
	x = random_vector (size1);
	b = c_matrix_dot_vector (1., a, x);
	c_vector_free (x);

	lu = c_matrix_alloc (a->size1, a->size2);
	c_matrix_memcpy (lu, a);

	x = c_vector_alloc (b->size);
	c_vector_memcpy (x, b);
	c_linalg_LU_solve (lu, x, NULL);
	c_matrix_free (lu);

	y = c_matrix_dot_vector (1., a, x);
	c_matrix_free (a);
	c_vector_free (x);

	/* b = -y + b */
	c_vector_axpy (-1., y, b);
	c_vector_free (y);

	nrm = c_vector_nrm (b);
	c_vector_free (b);

	return (nrm < 1.e-8);
}

bool
test_LU_svx (void)
{
	c_matrix		*a;
	c_vector		*x;
	c_vector		*b;

	c_matrix		*lu;
	c_vector_int	*p;
	c_vector		*y;
	double			nrm;

	a = random_matrix (size1, size1);
	x = random_vector (size1);
	b = c_matrix_dot_vector (1., a, x);
	c_vector_free (x);

	lu = c_matrix_alloc (a->size1, a->size2);
	c_matrix_memcpy (lu, a);
	c_linalg_LU_decomp (lu, &p);

	x = c_vector_alloc (b->size);
	c_vector_memcpy (x, b);
	c_linalg_LU_svx (lu, x, p);
	c_matrix_free (lu);
	c_vector_int_free (p);

	y = c_matrix_dot_vector (1., a, x);
	c_matrix_free (a);
	c_vector_free (x);

	/* b = -y + b */
	c_vector_axpy (-1., y, b);
	c_vector_free (y);

	nrm = c_vector_nrm (b);
	c_vector_free (b);

	return (nrm < 1.e-8);
}

bool
test_LU_invert (void)
{
	c_matrix		*a;
	c_matrix		*lu;
	c_matrix		*c;
	c_vector_int	*p;
	double			nrm;

	a = random_matrix (size1, size1);

	lu = c_matrix_alloc (a->size1, a->size2);
	c_matrix_memcpy (lu, a);
	c_linalg_LU_decomp (lu, &p);
	c_linalg_LU_invert (lu, p);
	c_vector_int_free (p);

	c = c_matrix_dot_matrix (1., lu, a);
	c_matrix_free (a);
	c_matrix_free (lu);

	nrm = c_matrix_nrm (c, '1');
	c_matrix_free (c);

	return (fabs (nrm - 1.) < 1.e-8);
}

bool
test_LU_1up (void)
{
	c_matrix		*a;
	c_vector		*t;
	c_vector		*s;

	c_matrix		*l;
	c_matrix		*u;
	c_vector_int	*p;
	double			nrm;

	a = random_matrix (size1, size2);
	{
		c_matrix		*lu = c_matrix_alloc (a->size1, a->size2);
		c_matrix_memcpy (lu, a);
		c_linalg_LU_decomp (lu, &p);
		c_linalg_LU_unpack (lu, &l, &u);
		c_matrix_free (lu);
	}

	s = random_vector (a->size1);
	t = random_vector (a->size2);
	{
		c_matrix	*s0 = c_matrix_view_array (s->size, 1, s->size, s->data);
		c_matrix	*t0 = c_matrix_view_array (t->size, 1, t->size, t->data);
		c_matrix	*c = c_matrix_dot_matrix_transpose (1., s0, t0);
		c_matrix_free (s0);
		c_matrix_free (t0);
		c_matrix_add (a, c);
		c_matrix_free (c);
	}

	c_linalg_LU_1up (l, u, p, s, t);
	c_vector_free (s);
	c_vector_free (t);

	{
		c_matrix	*a1 = c_matrix_dot_matrix (1., l, u);
		c_matrix_free (l);
		c_matrix_free (u);

		c_matrix_permute_rows (a, p);
		c_vector_int_free (p);
		c_matrix_sub (a, a1);
		c_matrix_free (a1);
	}
	nrm = c_matrix_nrm (a, '1');
	c_matrix_free (a);

	return nrm;
}
