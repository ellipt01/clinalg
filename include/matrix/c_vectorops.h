/*
 * c_vectorpos.h
 *
 *  Created on: 2014/04/10
 *      Author: utsugi
 */

#ifndef C_VECTORPOS_H_
#define C_VECTORPOS_H_

#ifdef __cplusplus
extern "C" {
#endif

void		c_vector_add_constant (c_vector *v, const double x);
double		c_vector_mean (const c_vector *v);
void		c_vector_sub (c_vector *y, const c_vector *x);
double		c_vector_asum (const c_vector *v);
int			c_vector_amax (const c_vector *v);
void		c_vector_scale (c_vector *v, double alpha);
double		c_vector_nrm (const c_vector *v);
void		c_vector_axpy (double alpha, const c_vector *x, c_vector *y);
double		c_vector_dot_vector (const c_vector *x, const c_vector *y);

#ifdef __cplusplus
}
#endif

#endif /* C_VECTORPOS_H_ */
