/*
 * test_chol.c
 *
 *  Created on: 2014/04/14
 *      Author: utsugi
 */

#include <clinalg.h>

#include "test_clinalg.h"

extern int		size1;

/* check |a - l' * l| < 1.e-8 */
bool
test_cholesky_decomp (void)
{
	c_matrix	*a;
	c_matrix	*c;
	c_matrix	*l;
	c_matrix	*b;
	double		nrm;

	{
		int			i;
		c_matrix	*a0 = random_matrix (size1, size1);
		a = c_matrix_transpose_dot_matrix (1., a0, a0);
		for (i = 0; i < size1; i++) c_matrix_set (a, i, i, c_matrix_get(a, i, i) + 0.1);
		c_matrix_free (a0);
	}

	/* c = chol(a) */
	c = c_matrix_alloc (a->size1, a->size2);
	c_matrix_memcpy (c, a);
	c_linalg_cholesky_decomp (c);
	l = c_matrix_alloc (c->size1, c->size2);
	c_matrix_set_zero (l);
	c_matrix_upper_triangular_memcpy (l, c);
	c_matrix_free (c);

	/* b = l' * l */
	b = c_matrix_transpose_dot_matrix (1., l, l);
	c_matrix_free (l);
	c_matrix_sub (a, b);
	c_matrix_free (b);

	nrm = c_matrix_nrm (a, '1');
	c_matrix_free (a);

	return (nrm < 1.e-8);
}

/* check |x - (l' * l)^-1 * y| < 1.e-8 */
bool
test_cholesky_svx (void)
{
	c_matrix	*a;
	c_vector	*x;
	c_vector	*y;

	c_matrix	*l;
	double		nrm;

	/* posdef symmetry matrix *a */
	{
		int			i;
		c_matrix	*a0 = random_matrix (size1, size1);
		a = c_matrix_transpose_dot_matrix (1., a0, a0);
		for (i = 0; i < size1; i++) c_matrix_set (a, i, i, c_matrix_get(a, i, i) + 0.1);
		c_matrix_free (a0);
	}
	/* vector *x */
	x = random_vector (size1);

	/* vector *y */
	y = c_matrix_dot_vector (1., a, x);

	/* cholesky_svx */
	c_linalg_cholesky_decomp (a);
	c_linalg_cholesky_svx (a, y);
	c_matrix_free (a);

	/* x = - y + x */
	c_vector_axpy (-1., y, x);
	c_vector_free (y);
	nrm = c_vector_nrm (x);
	c_vector_free (x);

	return (nrm < 1.e-8);
}

bool
test_cholesky_1up (void)
{
	c_matrix	*a;
	c_matrix	*c;
	c_matrix	*l;
	c_vector	*u;
	double		nrm;

	/* posdef symmetry matrix *a */
	{
		int			i;
		c_matrix	*a0 = random_matrix (size1, size1);
		a = c_matrix_transpose_dot_matrix (1., a0, a0);
		c_matrix_free (a0);
		for (i = 0; i < size1; i++) c_matrix_set (a, i, i, c_matrix_get(a, i, i) + 1.);
	}

	l = c_matrix_alloc (a->size1, a->size2);
	c_matrix_memcpy (l, a);

	u = random_vector (size1);
	{
		c_matrix	*ut = c_matrix_view_array (u->size, 1, u->size, u->data);
		c = c_matrix_dot_matrix_transpose (1., ut, ut);
		c_matrix_free (ut);
		c_matrix_add (a, c);
		c_matrix_free (c);
		c_linalg_cholesky_decomp (a);
	}

	c_linalg_cholesky_decomp (l);
	c_linalg_cholesky_1up (l, u);
	c_matrix_sub (a, l);
	c_matrix_free (l);

	nrm = c_matrix_nrm (a, '1');
	c_matrix_free (a);

	return (nrm < 1.e-8);
}

bool
test_cholesky_1down (void)
{
	int			info;
	c_matrix	*a;
	c_matrix	*c;
	c_matrix	*l;
	c_vector	*u;
	double		nrm;

	/* posdef symmetry matrix *a */
	{
		int			i;
		c_matrix	*a0 = random_matrix (size1, size1);
		a = c_matrix_transpose_dot_matrix (1., a0, a0);
		c_matrix_free (a0);
		for (i = 0; i < size1; i++) c_matrix_set (a, i, i, c_matrix_get(a, i, i) + 1.);
	}

	l = c_matrix_alloc (a->size1, a->size2);
	c_matrix_memcpy (l, a);
	u = random_vector (size1);
	c_vector_scale (u, 0.1);
	{
		c_matrix	*ut = c_matrix_view_array (u->size, 1, u->size, u->data);
		c = c_matrix_dot_matrix_transpose (1., ut, ut);
		c_matrix_free (ut);
		c_matrix_sub (a, c);
		c_matrix_free (c);
		c_linalg_cholesky_decomp (a);
	}

	c_linalg_cholesky_decomp (l);
	info = c_linalg_cholesky_1down (l, u);
	c_matrix_sub (a, l);
	c_matrix_free (l);

	nrm = c_matrix_nrm (a, '1');
	c_matrix_free (a);

	return (info == 0 && nrm < 1.e-8);
}

bool
test_cholesky_insert (void)
{
	int			index = size1 * rand () / RAND_MAX;
	c_matrix	*a;
	c_matrix	*b;
	c_matrix	*l;
	c_vector	*c;
	double		nrm;

	/* posdef symmetry matrix *a */
	{
		int			i;
		c_matrix	*a0 = random_matrix (size1, size1);
		a = c_matrix_transpose_dot_matrix (1., a0, a0);
		for (i = 0; i < size1; i++) c_matrix_set (a, i, i, c_matrix_get(a, i, i) + 0.1);
		c_matrix_free (a0);
	}

	l = c_matrix_alloc (size1 - 1, size1 - 1);
	c = c_vector_alloc (size1);
	{
		int			i, j, m, n;
		for (i = 0, m = 0; i < size1; i++) {
			c_vector_set (c, i, c_matrix_get (a, i, index));
			if (i == index) continue;
			for (j = 0, n = 0; j < size1; j++) {
				if (j == index) continue;
				c_matrix_set (l, m, n, c_matrix_get (a, i, j));
				n++;
			}
			m++;
		}
	}

	c_linalg_cholesky_decomp (l);
	c_linalg_cholesky_insert (l, index, c);
	c_vector_free (c);
	{
		int		i, j;
		for (i = 0; i < size1; i++) {
			for (j = 0; j < i; j++) c_matrix_set (l, i, j, 0.);
		}
	}
	b = c_matrix_transpose_dot_matrix (1., l, l);
	c_matrix_free (l);

	c_matrix_sub (a, b);
	c_matrix_free (b);

	nrm = c_matrix_nrm (a, '1');
	c_matrix_free (a);

	return (nrm < 1.e-8);
}

bool
test_cholesky_delete (void)
{
	int			index = size1 * rand () / RAND_MAX;
	c_matrix	*a;
	c_matrix	*l;
	c_matrix	*b;
	double		nrm;

	/* posdef symmetry matrix *a */
	{
		int			i;
		c_matrix	*a0 = random_matrix (size1, size1);
		l = c_matrix_transpose_dot_matrix (1., a0, a0);
		for (i = 0; i < size1; i++) c_matrix_set (l, i, i, c_matrix_get(l, i, i) + 0.1);
		c_matrix_free (a0);
	}

	a = c_matrix_alloc (size1 - 1, size1 - 1);
	{
		int			i, j, m, n;
		for (i = 0, m = 0; i < size1; i++) {
			if (i == index) continue;
			for (j = 0, n = 0; j < size1; j++) {
				if (j == index) continue;
				c_matrix_set (a, m, n, c_matrix_get (l, i, j));
				n++;
			}
			m++;
		}
	}

	c_linalg_cholesky_decomp (l);
	c_linalg_cholesky_delete (l, index);
	{
		int		i, j;
		for (i = 0; i < size1 - 1; i++) {
			for (j = 0; j < i; j++) c_matrix_set (l, i, j, 0.);
		}
	}
	b = c_matrix_transpose_dot_matrix (1., l, l);
	c_matrix_free (l);

	c_matrix_sub (a, b);
	c_matrix_free (b);

	nrm = c_matrix_nrm (a, '1');
	c_matrix_free (a);

	return (nrm < 1.e-8);
}
