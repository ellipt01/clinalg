/*
 * c_linalg_lu.h
 *
 *  Created on: 2014/04/10
 *      Author: utsugi
 */

#ifndef C_LINALG_LU_H_
#define C_LINALG_LU_H_

int			c_linalg_lapack_dgetrf (c_matrix *a, long **p);
int			c_linalg_lapack_dgetrs (char trans, c_matrix *lu, c_matrix *b, long *p);
int			c_linalg_lapack_dgetri (c_matrix *a, c_vector_int *p);
int			c_linalg_LU_decomp (c_matrix *a, long **p);
int			c_linalg_LU_solve (c_matrix *a, c_vector *b, long *p);
int			c_linalg_LU_invert (c_matrix *a);

#endif /* C_LINALG_LU_H_ */
