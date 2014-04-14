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
double			c_matrix_nrm (c_matrix *a, char norm);
c_matrix		*c_matrix_copy_upper_triangular (c_matrix *a);
c_matrix		*c_matrix_copy_lower_triangular (c_matrix *a);
c_vector		*c_matrix_dot_vector (double alpha, const c_matrix *a, const c_vector *v, double beta);
c_vector		*c_matrix_transpose_dot_vector (double alpha, const c_matrix *a, const c_vector *x, double beta);
c_matrix		*c_matrix_dot_matrix (double alpha, const c_matrix *a, const c_matrix *b, double beta);
c_matrix		*c_matrix_dot_matrix_transpose (double alpha, const c_matrix *a, const c_matrix *b, double beta);
c_matrix		*c_matrix_transpose_dot_matrix (double alpha, const c_matrix *a, const c_matrix *b, double beta);

#ifdef __cplusplus
}
#endif

#endif /* C_MATRIXOPS_H_ */
