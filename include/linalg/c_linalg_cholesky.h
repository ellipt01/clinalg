/*
 * c_linalg.h
 *
 *  Created on: 2014/04/03
 *      Author: utsugi
 */

#ifndef C_LINALG_CHOLESKY_H_
#define C_LINALG_CHOLESKY_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <c_matrix.h>

int		c_linalg_cholesky_decomp (c_matrix *a);
int		c_linalg_cholesky_solve (c_matrix *a, c_vector *b);
int		c_linalg_cholesky_svx (c_matrix *l, c_vector *b);
int		c_linalg_cholesky_invert (c_matrix *l);

void	c_linalg_cholesky_1up (c_matrix *l, c_vector *u);
int		c_linalg_cholesky_1down (c_matrix *l, c_vector *u);
int		c_linalg_cholesky_insert (c_matrix *l, const int index, c_vector *u);
void	c_linalg_cholesky_delete (c_matrix *l, const int index);

#ifdef __cplusplus
}
#endif

#endif /* C_LINALG_CHOLESKY_H_ */
