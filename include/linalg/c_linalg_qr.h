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

int		c_linalg_QR_decomp (c_matrix *a, c_vector_int **p, c_vector **tau);
int		c_linalg_QR_unpack (const c_matrix *qr, const c_vector *tau, c_matrix **q, c_matrix **r, bool econ);
int		c_linalg_QR_solve (c_matrix *a, c_vector *b);
int		c_linalg_lsQ_solve (double rcond, c_matrix *a, c_vector *b, c_vector_int **p, int *rank);
int		c_linalg_QR_Rsolve (c_matrix *r, c_vector *qty);
int		c_linalg_QR_RTsolve (c_matrix *r, c_vector *y);

void	c_linalg_QR_1up (c_matrix *q, c_matrix *r, const c_vector *u, const c_vector *v);
void	c_linalg_QR_colinsert (c_matrix *q, c_matrix *r, const int index, const c_vector *u);
void	c_linalg_QR_rowinsert (c_matrix *q, c_matrix *r, const int index, const c_vector *u);
void	c_linalg_QR_coldelete (c_matrix *q, c_matrix *r, const int index);
void	c_linalg_QR_rowdelete (c_matrix *q, c_matrix *r, const int index);

#ifdef __cplusplus
}
#endif

#endif /* C_LINALG_QR_H_ */
