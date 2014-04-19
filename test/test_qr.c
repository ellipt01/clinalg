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

	size1 = 50;
	size2 = 60;

	a = random_matrix (size1, size2);

	/* qr = qr(a) */
	qr = c_matrix_alloc (size1, size2);
	c_matrix_memcpy (qr, a);
	c_linalg_QR_decomp (qr, NULL, &tau);
	c_linalg_QR_unpack (qr, tau, &q, &r, false);
	c_vector_free (tau);
	c_matrix_free (qr);
/*
	fprintf (stdout, "a = [\n");
	c_matrix_fprintf2 (stdout, a, "%f");
	fprintf (stdout, "]\n");
	fprintf (stdout, "r = [\n");
	c_matrix_fprintf2 (stdout, r, "%f");
	fprintf (stdout, "]\n");
	fprintf (stdout, "q = [\n");
	c_matrix_fprintf2 (stdout, q, "%f");
	fprintf (stdout, "]\n");
*/
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
test_QR_decomp_econ (void)
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

	size1 = 60;
	size2 = 50;

	a = random_matrix (size1, size2);

	/* qr = qr(a) */
	qr = c_matrix_alloc (size1, size2);
	c_matrix_memcpy (qr, a);
	c_linalg_QR_decomp (qr, NULL, &tau);
	c_linalg_QR_unpack (qr, tau, &q, &r, true);
	c_vector_free (tau);
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
		c_matrix	*qr = c_matrix_alloc (a->size1, a->size2);
		c_matrix_memcpy (qr, a);
		c_linalg_QR_decomp (qr, NULL, &tau);
		c_linalg_QR_unpack (qr, tau, &q, &r, false);
		c_vector_free (tau);
		c_matrix_free (qr);

		/* x = q' * y */
		x = c_matrix_transpose_dot_vector (1., q, y, 0.);
		c_matrix_free (q);

		/* x = r^-1 * (q' * y) */
		c_linalg_QR_Rsolve (r, x);
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

bool
test_QR_1up ()
{
	double		nrm;
	size_t		size1 = 5;
	size_t		size2 = 6;
	c_matrix	*a;
	c_vector	*u;
	c_vector	*v;

	c_matrix	*q;
	c_matrix	*r;

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
		c_matrix	*c = c_matrix_dot_matrix_transpose (1., u0, v0, 0.);
		c_matrix_free (u0);
		c_matrix_free (v0);
		c_matrix_add (a, c);
		c_matrix_free (c);
	}

	c_linalg_QR_1up (q, r, u, v);
	c_vector_free (u);
	c_vector_free (v);

	{
		c_matrix	*a2 = c_matrix_dot_matrix (1., q, r, 0.);
		c_matrix_free (q);
		c_matrix_free (r);

		c_matrix_sub (a, a2);
		c_matrix_free (a2);
	}
	nrm = c_matrix_nrm (a, '1');
	c_matrix_free (a);

	return nrm;
}
