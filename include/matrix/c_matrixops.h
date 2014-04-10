/*
 * c_matrixops.h
 *
 *  Created on: 2014/04/10
 *      Author: utsugi
 */

#ifndef C_MATRIXOPS_H_
#define C_MATRIXOPS_H_

#ifdef __cplusplus
extern "C" {
#endif

/* c_matrixops.c */
c_vector		*c_matrix_dot_vector (double alpha, const c_matrix *a, const c_vector *v, double beta);
c_vector		*c_matrix_transpose_dot_vector (double alpha, const c_matrix *a, const c_vector *x, double beta);

#ifdef __cplusplus
}
#endif

#endif /* C_MATRIXOPS_H_ */
