/*
 * c_linalg_qr.h
 *
 *  Created on: 2014/04/11
 *      Author: utsugi
 */

#ifndef C_LINALG_QR_H_
#define C_LINALG_QR_H_

#ifdef __cplusplus
extern "C" {
#endif

int			c_linalg_lapack_dgeqrf (c_matrix *a, c_vector **tau);
int			c_linalg_lapack_dgeqp3 (c_matrix *a, c_vector **tau, c_vector_int **p);
int			c_linalg_lapack_dorgqr (c_matrix *qr, const c_vector *tau);
int			c_linalg_lapack_dgels (char trans, c_matrix *qr, c_matrix *b);
int			c_linalg_lapack_dgelsy (double rcond, c_matrix *qr, c_matrix *b, c_vector_int **p, int *rank);

int			c_linalg_QR_decomp (c_matrix *a, c_vector_int **p, c_vector **tau);
int			c_linalg_QR_unpack (const c_matrix *qr, const c_vector *tau, c_matrix **q, c_matrix **r, bool econ);
int			c_linalg_QR_solve (c_matrix *qr, c_vector *b);
int			c_linalg_lsQ_solve (double rcond, c_matrix *qr, c_vector *b, c_vector_int **p, int *rank);
void		c_linalg_QR_Rsolve (c_matrix *r, c_vector *qty);

void		c_linalg_QR_1up (c_matrix *q, c_matrix *r, const c_vector *u, const c_vector *v);
void		c_linalg_QR_colinsert (c_matrix *q, c_matrix *r, const size_t index, const c_vector *u);
void		c_linalg_QR_rowinsert (c_matrix *q, c_matrix *r, const size_t index, const c_vector *u);
void		c_linalg_QR_coldelete (c_matrix *q, c_matrix *r, const size_t index);
void		c_linalg_QR_rowdelete (c_matrix *q, c_matrix *r, const size_t index);

#ifdef __cplusplus
}
#endif

#endif /* C_LINALG_QR_H_ */
