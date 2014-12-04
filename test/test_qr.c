/*
 * test_qr.c
 *
 *  Created on: 2014/04/14
 *      Author: utsugi
 */

#include <clinalg.h>

#include "test_clinalg.h"

extern int		size1;
extern int		size2;

/* check |a - l' * l| < 1.e-8 */
bool
test_QR_decomp (void)
{
	c_matrix		*a;
	c_matrix		*qr;
	c_matrix		*q;
	c_matrix		*r;
	c_matrix		*b;
	c_vector		*tau;
	double			nrm;
	c_vector_int	*p;

	a = random_matrix (size1, size2);

	/* qr = qr(a) */
	qr = c_matrix_alloc (size1, size2);
	c_matrix_memcpy (qr, a);
	c_linalg_QR_decomp (qr, &p, &tau);
	c_linalg_QR_unpack (qr, tau, &q, &r, false);
	c_vector_free (tau);
	c_matrix_free (qr);

	/* b = q * r * p */
	b = c_matrix_dot_matrix (1., q, r);
	c_matrix_free (q);
	c_matrix_free (r);
	c_matrix_permute_cols (b, p);
	c_vector_int_free (p);

	c_matrix_sub (a, b);
	c_matrix_free (b);

	nrm = c_matrix_nrm (a, '1');
	c_matrix_free (a);

	return (nrm < 1.e-8);
}

bool
test_QR_decomp_econ (void)
{
	c_matrix	*a;
	c_matrix	*qr;
	c_matrix	*q;
	c_matrix	*r;
	c_matrix	*b;
	c_vector	*tau;
	double		nrm;

	a = random_matrix (size1, size2);

	/* qr = qr(a) */
	qr = c_matrix_alloc (size1, size2);
	c_matrix_memcpy (qr, a);
	c_linalg_QR_decomp (qr, NULL, &tau);
	c_linalg_QR_unpack (qr, tau, &q, &r, true);
	c_vector_free (tau);
	c_matrix_free (qr);

	/* b = q * r */
	b = c_matrix_dot_matrix (1., q, r);
	c_matrix_free (q);
	c_matrix_free (r);

	c_matrix_sub (a, b);
	c_matrix_free (b);

	nrm = c_matrix_nrm (a, '1');
	c_matrix_free (a);

	return (nrm < 1.e-8);
}

bool
test_QR_solve (void)
{
	c_matrix	*a;
	c_vector	*x;
	c_vector	*y;
	c_vector	*z;
	double		nrm;

	a = random_matrix (size1, size2);
	x = random_vector (size2);
	y = c_matrix_dot_vector (1., a, x);
	c_vector_free (x);

	{
		c_matrix	*tmp = c_matrix_alloc (a->size1, a->size2);
		c_matrix_memcpy (tmp, a);

		x = c_vector_alloc (y->size);
		c_vector_memcpy (x, y);

		c_linalg_QR_solve (tmp, x);
		c_matrix_free (tmp);
	}
	z = c_matrix_dot_vector (1., a, x);

	c_vector_axpy (-1., y, z);
	c_vector_free (y);

	nrm = c_vector_nrm (z);
	c_vector_free (z);

	return (nrm < 1.e-8);
}

bool
test_lsQ_solve (void)
{
	c_matrix	*a;
	c_vector	*x;
	c_vector	*y;
	c_vector	*z;
	double		nrm;

	a = random_matrix (size1, size2);
	x = random_vector (size2);
	y = c_matrix_dot_vector (1., a, x);
	c_vector_free (x);

	{
		c_matrix	*tmp = c_matrix_alloc (a->size1, a->size2);
		c_matrix_memcpy (tmp, a);

		x = c_vector_alloc (y->size);
		c_vector_memcpy (x, y);
		c_linalg_lsQ_solve (1.e-8, tmp, x, NULL, NULL);
		c_matrix_free (tmp);
	}
	z = c_matrix_dot_vector (1., a, x);

	c_vector_axpy (-1., y, z);
	c_vector_free (y);

	nrm = c_vector_nrm (z);
	c_vector_free (z);

	return (nrm < 1.e-8);
}

bool
test_QR_Rsolve (void)
{
	int 		_size1 = size2 + 10;
	c_matrix	*a;
	c_vector	*x;
	c_vector	*y;
	c_vector	*z;
	double		nrm;

	a = random_matrix (_size1, size2);
	x = random_vector (a->size2);
	y = c_matrix_dot_vector (1., a, x);
	c_vector_free (x);

	{
		c_vector	*tau;
		c_matrix	*q;
		c_matrix	*r;
		c_matrix	*qr = c_matrix_alloc (a->size1, a->size2);
		c_matrix_memcpy (qr, a);
		c_linalg_QR_decomp (qr, NULL, &tau);
		c_linalg_QR_unpack (qr, tau, &q, &r, true);
		c_vector_free (tau);
		c_matrix_free (qr);

		/* x = q' * y */
		x = c_matrix_transpose_dot_vector (1., q, y);
		c_matrix_free (q);

		/* x = r^-1 * (q' * y) */
		c_linalg_QR_Rsolve (r, x);
		c_matrix_free (r);
	}
	z = c_matrix_dot_vector (1., a, x);
	c_matrix_free (a);
	c_vector_free (x);

	c_vector_axpy (-1., y, z);
	c_vector_free (y);

	nrm = c_vector_nrm (z);
	c_vector_free (z);

	return (nrm < 1.e-8);
}

bool
test_QR_1up (void)
{
	c_matrix	*a;
	c_vector	*u;
	c_vector	*v;

	c_matrix	*q;
	c_matrix	*r;
	double		nrm;

	a = random_matrix (size1, size2);
	{
		c_vector	*tau;
		c_matrix	*qr = c_matrix_alloc (a->size1, a->size2);
		c_matrix_memcpy (qr, a);
		c_linalg_QR_decomp (qr, NULL, &tau);
		c_linalg_QR_unpack (qr, tau, &q, &r, false);
		c_vector_free (tau);
		c_matrix_free (qr);
	}

	u = random_vector (size1);
	v = random_vector (size2);
	{
		c_matrix	*u0 = c_matrix_view_array (u->size, 1, u->size, u->data);
		c_matrix	*v0 = c_matrix_view_array (v->size, 1, v->size, v->data);
		c_matrix	*c = c_matrix_dot_matrix_transpose (1., u0, v0);
		c_matrix_free (u0);
		c_matrix_free (v0);
		c_matrix_add (a, c);
		c_matrix_free (c);
	}

	c_linalg_QR_1up (q, r, u, v);
	c_vector_free (u);
	c_vector_free (v);

	{
		c_matrix	*a2 = c_matrix_dot_matrix (1., q, r);
		c_matrix_free (q);
		c_matrix_free (r);

		c_matrix_sub (a, a2);
		c_matrix_free (a2);
	}
	nrm = c_matrix_nrm (a, '1');
	c_matrix_free (a);

	return nrm;
}

bool
test_QR_colinsert (void)
{
	double		nrm;
	int			index = size2 * rand () / RAND_MAX;
	c_matrix	*a;
	c_vector	*u;

	c_matrix	*q;
	c_matrix	*r;

	c_matrix	*a1;
	c_matrix	*a2;

	a = random_matrix (size1, size2);
	u = random_vector (size1);

	/* QR decomposition */
	{
		c_vector	*tau;
		c_matrix	*qr = c_matrix_alloc (a->size1, a->size2);
		c_matrix_memcpy (qr, a);
		c_linalg_QR_decomp (qr, NULL, &tau);
		c_linalg_QR_unpack (qr, tau, &q, &r, true);
		c_vector_free (tau);
		c_matrix_free (qr);
	}
	c_linalg_QR_colinsert (q, r, index, u);
	a2 = c_matrix_dot_matrix (1., q, r);
	c_matrix_free (q);
	c_matrix_free (r);

	a1 = c_matrix_alloc (a->size1, a->size2 + 1);
	{
		int			j;
		c_vector	*col = c_vector_alloc (a->size1);

		for (j = 0; j < index; j++) {
			c_matrix_get_col (col, a, j);
			c_matrix_set_col (a1, j, col);
		}
		c_matrix_set_col (a1, index, u);
		c_vector_free (u);
		for (j = index; j < a->size2; j++) {
			c_matrix_get_col (col, a, j);
			c_matrix_set_col (a1, j + 1, col);
		}
		c_vector_free (col);
		c_matrix_free (a);
	}

	c_matrix_sub (a1, a2);
	c_matrix_free (a2);

	nrm = c_matrix_nrm (a1, '1');
	c_matrix_free (a1);

	return nrm;
}

bool
test_QR_rowinsert (void)
{
	double		nrm;
	int			index = size1 * rand () / RAND_MAX;
	c_matrix	*a;
	c_vector	*u;

	c_matrix	*q;
	c_matrix	*r;

	c_matrix	*a1;
	c_matrix	*a2;

	a = random_matrix (size1, size2);
	u = random_vector (size2);

	/* QR decomposition */
	{
		c_vector	*tau;
		c_matrix	*qr = c_matrix_alloc (a->size1, a->size2);
		c_matrix_memcpy (qr, a);
		c_linalg_QR_decomp (qr, NULL, &tau);
		c_linalg_QR_unpack (qr, tau, &q, &r, false);
		c_vector_free (tau);
		c_matrix_free (qr);
	}
	c_linalg_QR_rowinsert (q, r, index, u);
	a2 = c_matrix_dot_matrix (1., q, r);
	c_matrix_free (q);
	c_matrix_free (r);

	a1 = c_matrix_alloc (a->size1 + 1, a->size2);
	{
		int			i;
		c_vector	*row = c_vector_alloc (a->size2);

		for (i = 0; i < index; i++) {
			c_matrix_get_row (row, a, i);
			c_matrix_set_row (a1, i, row);
		}
		c_matrix_set_row (a1, index, u);
		c_vector_free (u);
		for (i = index; i < a->size1; i++) {
			c_matrix_get_row (row, a, i);
			c_matrix_set_row (a1, i + 1, row);
		}
		c_vector_free (row);
		c_matrix_free (a);
	}

	c_matrix_sub (a1, a2);
	c_matrix_free (a2);

	nrm = c_matrix_nrm (a1, '1');
	c_matrix_free (a1);

	return nrm;
}

bool
test_QR_coldelete (void)
{
	double		nrm;
	int			index = size2 * rand () / RAND_MAX;
	c_matrix	*a;
	c_matrix	*q;
	c_matrix	*r;

	c_matrix	*a1;
	c_matrix	*a2;

	a = random_matrix (size1, size2);

	/* QR decomposition */
	{
		c_vector	*tau;
		c_matrix	*qr = c_matrix_alloc (a->size1, a->size2);
		c_matrix_memcpy (qr, a);
		c_linalg_QR_decomp (qr, NULL, &tau);
		c_linalg_QR_unpack (qr, tau, &q, &r, true);
		c_vector_free (tau);
		c_matrix_free (qr);
	}
	c_linalg_QR_coldelete (q, r, index);
	a2 = c_matrix_dot_matrix (1., q, r);
	c_matrix_free (q);
	c_matrix_free (r);

	a1 = c_matrix_alloc (a->size1, a->size2 - 1);
	{
		int			j;
		c_vector	*col = c_vector_alloc (a->size1);

		for (j = 0; j < index; j++) {
			c_matrix_get_col (col, a, j);
			c_matrix_set_col (a1, j, col);
		}
		for (j = index + 1; j < a->size2; j++) {
			c_matrix_get_col (col, a, j);
			c_matrix_set_col (a1, j - 1, col);
		}
		c_vector_free (col);
		c_matrix_free (a);
	}

	c_matrix_sub (a1, a2);
	c_matrix_free (a2);

	nrm = c_matrix_nrm (a1, '1');
	c_matrix_free (a1);

	return nrm;
}

bool
test_QR_rowdelete (void)
{
	double		nrm;
	int			index = size1 * rand () / RAND_MAX;
	c_matrix	*a;

	c_matrix	*q;
	c_matrix	*r;

	c_matrix	*a1;
	c_matrix	*a2;

	a = random_matrix (size1, size2);

	/* QR decomposition */
	{
		c_vector	*tau;
		c_matrix	*qr = c_matrix_alloc (a->size1, a->size2);
		c_matrix_memcpy (qr, a);
		c_linalg_QR_decomp (qr, NULL, &tau);
		c_linalg_QR_unpack (qr, tau, &q, &r, false);
		c_vector_free (tau);
		c_matrix_free (qr);
	}
	c_linalg_QR_rowdelete (q, r, index);
	a2 = c_matrix_dot_matrix (1., q, r);
	c_matrix_free (q);
	c_matrix_free (r);

	a1 = c_matrix_alloc (a->size1 - 1, a->size2);
	{
		int			i;
		c_vector	*row = c_vector_alloc (a->size2);

		for (i = 0; i < index; i++) {
			c_matrix_get_row (row, a, i);
			c_matrix_set_row (a1, i, row);
		}
		for (i = index + 1; i < a->size1; i++) {
			c_matrix_get_row (row, a, i);
			c_matrix_set_row (a1, i - 1, row);
		}
		c_vector_free (row);
		c_matrix_free (a);
	}

	c_matrix_sub (a1, a2);
	c_matrix_free (a2);

	nrm = c_matrix_nrm (a1, '1');
	c_matrix_free (a1);

	return nrm;
}
