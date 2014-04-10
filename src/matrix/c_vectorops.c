/*
 * c_vectorops.c
 *
 *  Created on: 2014/04/03
 *      Author: utsugi
 */

#include <c_matrix.h>

/* c_matrix.c */
extern void	c_error (const char * function_name, const char *error_msg);

/* blas */
extern long	idamax_ (long *n, double *x, long *incx);
extern double	dasum_ (long *n, double *x, long *incx);
extern double	dnrm2_ (long *n, double *x, long *incx);
extern double	ddot_ (long *n, double *x, long *incx, double *y, long *incy);
extern void	dscal_ (long *n, double *alpha, double *x, long *incx);
extern void	daxpy_ (long *n, double *alpha, double *x, long *incx, double *y, long *incy);

void
c_vector_add_constant (c_vector *v, const double x)
{
	int		i;
	if (c_vector_is_empty (v)) c_error ("c_vector_add_constant", "vector is empty.");
	for (i = 0; i < v->size; i++) v->data[i * v->stride] += x;
	return;
}

double
c_vector_mean (const c_vector *v)
{
	int	i;
	double	sum;
	if (c_vector_is_empty (v)) c_error ("c_vector_mean", "vector is empty.");
	/* x = sum x / N */
	for (i = 0, sum = 0.0; i < v->size; i++) sum += v->data[i * v->stride];
	return sum / (double) v->size;
}

void
c_vector_sub (c_vector *y, const c_vector *x)
{
	long	n;
	long	incx;
	long	incy;
	double	alpha = - 1.;
	if (c_vector_is_empty (y)) c_error ("c_vector_sub", "first vector is empty.");
	if (c_vector_is_empty (x)) c_error ("c_vector_sub", "second vector is empty.");
	if (x->size != y->size) c_error ("c_vector_sub", "vector size does not match.");
	/* y = y - x */
	n = (long) x->size;
	incx = (long) x->stride;
	incy = (long) y->stride;
	daxpy_ (&n, &alpha, x->data, &incx, y->data, &incy);
	return;
}

double
c_vector_asum (const c_vector *v)
{
	long	n;
	long	incv;
	if (c_vector_is_empty (v)) c_error ("c_vector_asum", "vector is empty.");
	/* x = sum |x| */
	n = (long) v->size;
	incv = (long) v->stride;
	return dasum_ (&n, v->data, &incv);
}

int
c_vector_amax (const c_vector *v)
{
	long	n;
	long	incv;
	if (c_vector_is_empty (v)) c_error ("c_vector_amax", "vector is empty.");

	n = (long) v->size;
	incv = (long) v->stride;
	return (int) idamax_ (&n, v->data, &incv) - 1;
}

void
c_vector_scale (c_vector *v, double alpha)
{
	long	n;
	long	incv;
	if (c_vector_is_empty (v)) c_error ("c_vector_scale", "vector is empty.");
	/* x = alpha * x */
	n = (long) v->size;
	incv = (long) v->stride;
	dscal_ (&n, &alpha, v->data, &incv);
	return;
}

double
c_vector_nrm (const c_vector *v)
{
	long	n;
	long	incv;
	if (c_vector_is_empty (v)) c_error ("c_vector_nrm", "vector is empty.");
	n = (long) v->size;
	incv = (long) v->stride;
	return dnrm2_ (&n, v->data, &incv);
}

void
c_vector_axpy (double alpha, const c_vector *x, c_vector *y)
{
	long	n;
	long	incx;
	long	incy;
	if (c_vector_is_empty (x)) c_error ("c_vector_axpy", "first vector is empty.");
	if (c_vector_is_empty (y)) c_error ("c_vector_axpy", "second vector is empty.");
	if (x->size != y->size) c_error ("c_vector_axpy", "vector size does not match.");
	/* y = y + alpha * x */
	n = (long) x->size;
	incx = (long) x->stride;
	incy = (long) y->stride;
	daxpy_ (&n, &alpha, x->data, &incx, y->data, &incy);
	return;
}

double
c_vector_dot_vector (const c_vector *x, const c_vector *y)
{
	long	n;
	long	incx;
	long	incy;
	if (c_vector_is_empty (x)) c_error ("c_vector_dot_vector", "first vector is empty.");
	if (c_vector_is_empty (y)) c_error ("c_vector_dot_vector", "second vector is empty.");
	if (x->size != y->size) c_error ("c_vector_dot_vector", "vector size does not match.");
	n = (long) x->size;
	incx = (long) x->stride;
	incy = (long) y->stride;
	return ddot_ (&n, x->data, &incx, y->data, &incy);
}