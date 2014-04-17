/*
 * c_linalg_utils.h
 *
 *  Created on: 2014/04/17
 *      Author: utsugi
 */

#ifndef C_LINALG_UTILS_H_
#define C_LINALG_UTILS_H_

c_matrix	*c_linalg_permutation_matrix_row (const size_t size1, const size_t size2, const c_vector_int *p);
c_matrix	*c_linalg_permutation_matrix_col (const size_t size1, const size_t size2, const c_vector_int *p);


#endif /* C_LINALG_UTILS_H_ */
