/*
 * c_linalg_lapack.h
 *
 *  Created on: 2014/05/13
 *      Author: utsugi
 */

#ifndef C_LINALG_LAPACK_H_
#define C_LINALG_LAPACK_H_

#ifdef __cplusplus
extern "C" {
#endif

int		c_linalg_lapack_dtrtrs (char uplo, char trans, char diag, c_matrix *a, c_matrix *b);

/* LU decomposition */
int		c_linalg_lapack_dgetrf (c_matrix *a, c_vector_int **p);
int		c_linalg_lapack_dgetrs (char trans, c_matrix *lu, c_matrix *b, c_vector_int *p);
int		c_linalg_lapack_dgesv (c_matrix *a, c_matrix *b, c_vector_int **p);
int		c_linalg_lapack_dgetri (c_matrix *lu, c_vector_int *p);

/* Cholesky decomposition */
int		c_linalg_lapack_dpotrf (char uplo, c_matrix *a);
int		c_linalg_lapack_dpotrs (char uplo, c_matrix *l, c_matrix *b);
int		c_linalg_lapack_dpotri (char uplo, c_matrix *l);

/* QR decomposition */
int			c_linalg_lapack_dgeqrf (c_matrix *a, c_vector **tau);
int			c_linalg_lapack_dgeqp3 (c_matrix *a, c_vector **tau, c_vector_int **p);
int			c_linalg_lapack_dorgqr (c_matrix *qr, const c_vector *tau);
int			c_linalg_lapack_dgels (char trans, c_matrix *a, c_matrix *b);
int			c_linalg_lapack_dgelsy (double rcond, c_matrix *a, c_matrix *b, c_vector_int **p, int *rank);

/* SV decomposition */
int		c_linalg_lapack_dgesvd (char jobu, char jobvt, c_matrix *a, c_matrix **u, c_matrix **vt, c_vector **s);
int		c_linalg_lapack_dgesdd (char jobz, c_matrix *a, c_matrix **u, c_matrix **vt, c_vector **s);
int		c_linalg_lapack_dgelss (double rcond, c_matrix *a, c_matrix *b, c_vector **s, int *rank);
int		c_linalg_lapack_dgelsd (double rcond, c_matrix *a, c_matrix *b, c_vector **s, int *rank);

#ifdef __cplusplus
}
#endif

#endif /* C_LINALG_LAPACK_H_ */
