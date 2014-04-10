/*
 * c_vector_double.h
 *
 *  Created on: 2014/04/03
 *      Author: utsugi
 */

#ifndef C_VECTOR_DOUBLE_H_
#define C_VECTOR_DOUBLE_H_

#ifdef __cplusplus
extern "C" {
#endif

typedef struct s_c_vector	c_vector;

struct s_c_vector {
	size_t		size;
	size_t		stride;
	size_t		tsize;
	bool		owner;
	double		*data;
};

c_vector		*c_vector_alloc (const size_t size);
c_vector		*c_vector_view_array (const size_t size, const size_t stride, double *data);
bool			c_vector_is_empty (const c_vector *a);
void			c_vector_free (c_vector *v);
void			c_vector_set (c_vector *v, const int i, double val);
double			c_vector_get (const c_vector *v, const int i);
void			c_vector_memcpy (c_vector *dest, const c_vector *src);
void			c_vector_set_zero (c_vector *v);
void			c_vector_fprintf (FILE *stream, const c_vector *v, const char *format);

/* c_matrixops.c */
void			c_vector_add_constant (c_vector *v, const double x);
double			c_vector_mean (const c_vector *v);
void			c_vector_sub (c_vector *y, const c_vector *x);
double			c_vector_asum (const c_vector *v);
int				c_vector_amax (const c_vector *v);
void			c_vector_scale (c_vector *v, const double alpha);
double			c_vector_nrm (const c_vector *v);
void			c_vector_axpy (double alpha, const c_vector *x, c_vector *y);

#ifdef __cplusplus
}
#endif

#endif /* C_VECTOR_DOUBLE_H_ */
