/*
 * test_qr.c
 *
 *  Created on: 2014/04/14
 *      Author: utsugi
 */

#include <c_linalg.h>

#include "test_clinalg.h"

/* check |a - l' * l| < 1.e-8 */
bool
test_QR_decomp (void)
{
	size_t		size1;
	size_t		size2;
	c_matrix	*a;
	c_matrix	*qr;
	c_matrix	*q;
	c_matrix	*r;
	c_matrix	*b;
	c_vector	*tau;
	double		nrm;

	size1 = 5;
	size2 = 6;

	a = random_matrix (size1, size2);

	/* qr = qr(a) */
	qr = c_matrix_alloc (size1, size2);
	c_matrix_memcpy (qr, a);
	c_linalg_QR_decomp (qr, NULL, &tau);

	/* q */
	q = c_matrix_alloc (qr->size1, qr->size2);
	c_matrix_memcpy (q, qr);
	c_linalg_QR_unpack (q, tau);
	c_vector_free (tau);

	r = c_matrix_alloc (qr->size1, qr->size2);
	c_matrix_set_zero (r);
	c_matrix_upper_triangular_memcpy (r, qr);
	c_matrix_free (qr);

	/* b = q * r */
	b = c_matrix_dot_matrix (1., q, r, 0.);
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
	size_t		size1;
	size_t		size2;
	c_matrix	*a;
	c_vector	*x;
	c_vector	*y;
	c_vector	*z;
	double		nrm;

	size1 = 50;
	size2 = 60;

	a = random_matrix (size1, size2);
	x = random_vector (size2);
	y = c_matrix_dot_vector (1., a, x, 0.);
	c_vector_free (x);

	{
		c_matrix	*tmp = c_matrix_alloc (a->size1, a->size2);
		c_matrix_memcpy (tmp, a);

		x = c_vector_alloc (y->size);
		c_vector_memcpy (x, y);

		c_linalg_QR_solve (tmp, x);
		c_matrix_free (tmp);
	}
	z = c_matrix_dot_vector (1., a, x, 0.);

	c_vector_sub (z, y);
	c_vector_free (y);

	nrm = c_vector_nrm (z);
	c_vector_free (z);

	return (nrm < 1.e-8);
}

bool
test_lsQ_solve (void)
{
	size_t		size1;
	size_t		size2;
	c_matrix	*a;
	c_vector	*x;
	c_vector	*y;
	c_vector	*z;
	double		nrm;

	size1 = 5;
	size2 = 6;

	a = random_matrix (size1, size2);
	x = random_vector (size2);
	y = c_matrix_dot_vector (1., a, x, 0.);
	c_vector_free (x);

	{
		c_matrix	*tmp = c_matrix_alloc (a->size1, a->size2);
		c_matrix_memcpy (tmp, a);

		x = c_vector_alloc (y->size);
		c_vector_memcpy (x, y);
		c_linalg_lsQ_solve (1.e-8, tmp, x, NULL, NULL);
		c_matrix_free (tmp);
	}
	z = c_matrix_dot_vector (1., a, x, 0.);

	c_vector_sub (z, y);
	c_vector_free (y);

	nrm = c_vector_nrm (z);
	c_vector_free (z);

	return (nrm < 1.e-8);
}

bool
test_QR_Rsolve (void)
{
	size_t		size1;
	size_t		size2;
	c_matrix	*a;
	c_vector	*x;
	c_vector	*y;
	c_vector	*z;
	double		nrm;

	size1 = 60;
	size2 = 50;

	a = random_matrix (size1, size2);
	x = random_vector (size2);
	y = c_matrix_dot_vector (1., a, x, 0.);
	c_vector_free (x);

	{
		c_vector	*tau;
		c_matrix	*q;
		c_matrix	*r;
		c_matrix	*tmp = c_matrix_alloc (a->size1, a->size2);
		c_matrix_memcpy (tmp, a);

		c_linalg_QR_decomp (tmp, NULL, &tau);

		q = c_matrix_alloc (tmp->size1, tmp->size2);
		c_matrix_memcpy (q, tmp);
		c_linalg_QR_unpack (q, tau);
		c_vector_free (tau);

		/* x = q' * y */
		x = c_matrix_transpose_dot_vector (1., q, y, 0.);
		c_matrix_free (q);

		r = c_matrix_submatrix (0, 0, size2, size2, tmp);
		c_linalg_QR_Rsolve (r, x);
		c_matrix_free (tmp);
		c_matrix_free (r);
	}
	z = c_matrix_dot_vector (1., a, x, 0.);
	c_matrix_free (a);
	c_vector_free (x);

	c_vector_sub (z, y);
	c_vector_free (y);

	nrm = c_vector_nrm (z);
	c_vector_free (z);

	return (nrm < 1.e-8);
}
