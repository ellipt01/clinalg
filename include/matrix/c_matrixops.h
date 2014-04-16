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
void			c_matrix_swap_rows (const size_t i, const size_t j, c_matrix *a);
void			c_matrix_swap_cols (const size_t i, const size_t j, c_matrix *a);
double			c_matrix_nrm (c_matrix *a, char norm);
void			c_matrix_upper_triangular_memcpy (c_matrix *tr, const c_matrix *a);
void			c_matrix_lower_triangular_memcpy (c_matrix *tr, const c_matrix *a);
c_vector		*c_matrix_get_diagonal (const c_matrix *a);
c_matrix		*c_matrix_identity (const size_t size1, const size_t size2);
void			c_matrix_set_diagonal (const c_vector *d, c_matrix *a);
c_vector		*c_matrix_dot_vector (double alpha, const c_matrix *a, const c_vector *v, double beta);
c_vector		*c_matrix_transpose_dot_vector (double alpha, const c_matrix *a, const c_vector *x, double beta);
c_matrix		*c_matrix_dot_matrix (double alpha, const c_matrix *a, const c_matrix *b, double beta);
c_matrix		*c_matrix_dot_matrix_transpose (double alpha, const c_matrix *a, const c_matrix *b, double beta);
c_matrix		*c_matrix_transpose_dot_matrix (double alpha, const c_matrix *a, const c_matrix *b, double beta);

#ifdef __cplusplus
}
#endif

#endif /* C_MATRIXOPS_H_ */
