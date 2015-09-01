/*
 * c_linalg_eigen.h
 *
 *  Created on: 2015/09/01
 *      Author: utsugi
 */

#ifndef LINALG_C_LINALG_EIGEN_H_
#define LINALG_C_LINALG_EIGEN_H_

int		c_linalg_eigen_symm (c_matrix *a, c_vector **e);
int		c_linalg_eigen_symm_d (c_matrix *a, c_vector **e);
int		c_linalg_eigen_general (c_matrix *a, c_vector **er, c_vector **ei, c_matrix **vl, c_matrix **vr);
int		c_linalg_schur_decomp (c_matrix *a, c_vector **er, c_vector **ei, c_matrix **t);

#endif /* LINALG_C_LINALG_EIGEN_H_ */
