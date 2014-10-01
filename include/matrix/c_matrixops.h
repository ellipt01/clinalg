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
void			c_matrix_add (c_matrix *y, const c_matrix *x);
void			c_matrix_sub (c_matrix *y, const c_matrix *x);
void			c_matrix_axpy (double alpha, const c_matrix *x, c_matrix *y);
double			c_matrix_nrm (c_matrix *a, char norm);

void			c_matrix_swap_rows (const size_t i, const size_t j, c_matrix *a);
void			c_matrix_swap_cols (const size_t i, const size_t j, c_matrix *a);

void			c_matrix_permute_rows (c_matrix *a, const c_vector_int *p);
void			c_matrix_permute_cols (c_matrix *a, const c_vector_int *p);

void			c_matrix_add_rowcols (c_matrix *a, const size_t dm, const size_t dn);
void			c_matrix_remove_rowcols (c_matrix *a, const size_t dm, const size_t dn);

void			c_matrix_merge_row (c_matrix *a, const c_matrix *b);
void			c_matrix_merge_col (c_matrix *a, const c_matrix *b);

c_matrix		*c_matrix_identity (const size_t size1, const size_t size2);

c_vector		*c_matrix_dot_vector (double alpha, const c_matrix *a, const c_vector *v, double beta);
c_vector		*c_matrix_transpose_dot_vector (double alpha, const c_matrix *a, const c_vector *x, double beta);
c_vector		*c_matrix_symm_upper_dot_vector (double alpha, const c_matrix *a, const c_vector *x, double beta);
c_vector		*c_matrix_symm_lower_dot_vector (double alpha, const c_matrix *a, const c_vector *x, double beta);

c_matrix		*c_matrix_dot_matrix (double alpha, const c_matrix *a, const c_matrix *b, double beta);
c_matrix		*c_matrix_dot_matrix_transpose (double alpha, const c_matrix *a, const c_matrix *b, double beta);
c_matrix		*c_matrix_transpose_dot_matrix (double alpha, const c_matrix *a, const c_matrix *b, double beta);
c_matrix		*c_matrix_transpose_dot_matrix_transpose (double alpha, const c_matrix *a, const c_matrix *b, double beta);
c_matrix		*c_matrix_symm_upper_dot_matrix (double alpha, const c_matrix *a, const c_matrix *b, double beta);
c_matrix		*c_matrix_symm_lower_dot_matrix (double alpha, const c_matrix *a, const c_matrix *b, double beta);
c_matrix		*c_matrix_dot_matrix_symm_upper (double alpha, const c_matrix *a, const c_matrix *b, double beta);
c_matrix		*c_matrix_dot_matrix_symm_lower (double alpha, const c_matrix *a, const c_matrix *b, double beta);

#ifdef __cplusplus
}
#endif

#endif /* C_MATRIXOPS_H_ */
