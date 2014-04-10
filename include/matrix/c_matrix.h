/*
 * c_matrix.h
 *
 *  Created on: 2014/04/10
 *      Author: utsugi
 */

#ifndef C_MATRIX_H_
#define C_MATRIX_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include <c_matrix_double.h>
#include <c_matrixops.h>
#include <c_vector_double.h>
#include <c_vectorops.h>
#include <c_vector_int.h>

#define GET_INDEX_OF_VECTOR(v, i) (i * v->stride)
#define GET_INDEX_OF_MATRIX(a, i, j) (i + j * a->lda)

#ifdef __cplusplus
}
#endif

#endif /* C_MATRIX_H_ */
