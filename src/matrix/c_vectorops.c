/*
 * c_vectorops.c
 *
 *  Created on: 2014/04/03
 *      Author: utsugi
 */

#include <c_vector.h>
#include "../../include/clinalg_utils.h"
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

/* sum_i |x(i)| */
double
c_vector_asum (const c_vector *x)
{
	if (c_vector_is_empty (x)) c_error ("c_vector_asum", "vector is empty.");
	/* x = sum |x| */
	return F77CALL (dasum) (&x->size, x->data, &x->stride);
}

/* max_i |x(i)| */
int
c_vector_amax (const c_vector *x)
{
	if (c_vector_is_empty (x)) c_error ("c_vector_amax", "vector is empty.");
	return (int) F77CALL (idamax) (&x->size, x->data, &x->stride) - 1;
}

/* x = alpha * x */
void
c_vector_scale (c_vector *x, const double alpha)
{
	if (c_vector_is_empty (x)) c_error ("c_vector_scale", "vector is empty.");
	F77CALL (dscal) (&x->size, &alpha, x->data, &x->stride);
	return;
}

/* nrm (x) */
double
c_vector_nrm (const c_vector *x)
{
	if (c_vector_is_empty (x)) c_error ("c_vector_nrm", "vector is empty.");
	return F77CALL (dnrm2) (&x->size, x->data, &x->stride);
}

/* y = a * x + y */
void
c_vector_axpy (const double alpha, const c_vector *x, c_vector *y)
{
	if (c_vector_is_empty (x)) c_error ("c_vector_axpy", "first vector is empty.");
	if (c_vector_is_empty (y)) c_error ("c_vector_axpy", "second vector is empty.");
	if (x->size != y->size) c_error ("c_vector_axpy", "vector size does not match.");
	F77CALL (daxpy) (&x->size, &alpha, x->data, &x->stride, y->data, &y->stride);
	return;
}

/* x' * y */
double
c_vector_dot_vector (const c_vector *x, const c_vector *y)
{
	if (c_vector_is_empty (x)) c_error ("c_vector_dot_vector", "first vector is empty.");
	if (c_vector_is_empty (y)) c_error ("c_vector_dot_vector", "second vector is empty.");
	if (x->size != y->size) c_error ("c_vector_dot_vector", "vector size does not match.");
	return F77CALL (ddot) (&x->size, x->data, &x->stride, y->data, &y->stride);
}
