/*
 * c_vectorops.c
 *
 *  Created on: 2014/04/03
 *      Author: utsugi
 */

#include <c_vector.h>
#include <c_linalg_utils.h>

#include "private.h"

/* x = x + c */
void
c_vector_add_constant (c_vector *x, const double c)
{
	int		i;
	if (c_vector_is_empty (x)) c_error ("c_vector_add_constant", "vector is empty.");
	for (i = 0; i < x->size; i++) x->data[i * x->stride] += c;
	return;
}

/* sum_i x(i) */
double
c_vector_sum (const c_vector *x)
{
	int		i;
	double	sum;
	if (c_vector_is_empty (x)) c_error ("c_vector_mean", "vector is empty.");
	/* x = sum x / N */
	for (i = 0, sum = 0.0; i < x->size; i++) sum += x->data[i * x->stride];
	return sum;
}

/* sum_i x(i) / x->size */
double
c_vector_mean (const c_vector *x)
{
	double	sum;
	if (c_vector_is_empty (x)) c_error ("c_vector_mean", "vector is empty.");
	/* x = sum x / N */
	sum = c_vector_sum (x);
	return sum / (double) x->size;
}

/* y = y - x */
void
c_vector_sub (c_vector *y, const c_vector *x)
{
	int		n;
	int		incx;
	int		incy;
	double	alpha = - 1.;
	if (c_vector_is_empty (y)) c_error ("c_vector_sub", "first vector is empty.");
	if (c_vector_is_empty (x)) c_error ("c_vector_sub", "second vector is empty.");
	if (x->size != y->size) c_error ("c_vector_sub", "vector size does not match.");
	/* y = y - x */
	n = (int) x->size;
	incx = (int) x->stride;
	incy = (int) y->stride;
	daxpy_ (&n, &alpha, x->data, &incx, y->data, &incy);
	return;
}

/* sum_i |x(i)| */
double
c_vector_asum (const c_vector *x)
{
	int		n;
	int		incv;
	if (c_vector_is_empty (x)) c_error ("c_vector_asum", "vector is empty.");
	/* x = sum |x| */
	n = (int) x->size;
	incv = (int) x->stride;
	return dasum_ (&n, x->data, &incv);
}

/* max_i |x(i)| */
int
c_vector_amax (const c_vector *x)
{
	int		n;
	int		incv;
	if (c_vector_is_empty (x)) c_error ("c_vector_amax", "vector is empty.");

	n = (int) x->size;
	incv = (int) x->stride;
	return (int) idamax_ (&n, x->data, &incv) - 1;
}

/* x = alpha * x */
void
c_vector_scale (c_vector *x, const double alpha)
{
	int		n;
	int		incv;
	if (c_vector_is_empty (x)) c_error ("c_vector_scale", "vector is empty.");
	/* x = alpha * x */
	n = (int) x->size;
	incv = (int) x->stride;
	dscal_ (&n, &alpha, x->data, &incv);
	return;
}

/* sqrt (x' * x) */
double
c_vector_nrm (const c_vector *x)
{
	int		n;
	int		incv;
	if (c_vector_is_empty (x)) c_error ("c_vector_nrm", "vector is empty.");
	n = (int) x->size;
	incv = (int) x->stride;
	return dnrm2_ (&n, x->data, &incv);
}

/* y = a * x + y */
void
c_vector_axpy (const double alpha, const c_vector *x, c_vector *y)
{
	int		n;
	int		incx;
	int		incy;
	if (c_vector_is_empty (x)) c_error ("c_vector_axpy", "first vector is empty.");
	if (c_vector_is_empty (y)) c_error ("c_vector_axpy", "second vector is empty.");
	if (x->size != y->size) c_error ("c_vector_axpy", "vector size does not match.");
	/* y = y + alpha * x */
	n = (int) x->size;
	incx = (int) x->stride;
	incy = (int) y->stride;
	daxpy_ (&n, &alpha, x->data, &incx, y->data, &incy);
	return;
}

/* x' * y */
double
c_vector_dot_vector (const c_vector *x, const c_vector *y)
{
	int		n;
	int		incx;
	int		incy;
	if (c_vector_is_empty (x)) c_error ("c_vector_dot_vector", "first vector is empty.");
	if (c_vector_is_empty (y)) c_error ("c_vector_dot_vector", "second vector is empty.");
	if (x->size != y->size) c_error ("c_vector_dot_vector", "vector size does not match.");
	n = (int) x->size;
	incx = (int) x->stride;
	incy = (int) y->stride;
	return ddot_ (&n, x->data, &incx, y->data, &incy);
}
