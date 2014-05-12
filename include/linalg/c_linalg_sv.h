/*
 * c_linalg_sv.h
 *
 *  Created on: 2014/04/19
 *      Author: utsugi
 */

#ifndef C_LINALG_SV_H_
#define C_LINALG_SV_H_

#ifdef __cplusplus
extern "C" {
#endif

int		c_linalg_SV_decomp (c_matrix *a, c_matrix **u, c_matrix **vt, c_vector **s, bool econ);
int		c_linalg_SV_solve (double rcond, c_matrix *a, c_vector *b, c_vector **s, int *rank);
int		c_linalg_SV_lsd_solve (double rcond, c_matrix *a, c_vector *b, c_vector **s, int *rank);

#ifdef __cplusplus
}
#endif

#endif /* C_LINALG_SV_H_ */
