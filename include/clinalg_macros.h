/*
 * clinalg_macros.h
 *
 *  Created on: 2014/05/08
 *      Author: utsugi
 */

#ifndef CLINALG_MACROS_H_
#define CLINALG_MACROS_H_

#ifndef C_MAX
#define C_MAX(a, b)	((a) >= (b)) ? (a) : (b)
#endif

#ifndef C_MIN
#define C_MIN(a, b)	((a) >= (b)) ? (b) : (a)
#endif

#ifndef INDEX_OF_VECTOR
#define INDEX_OF_VECTOR(v, i) ((i) * (v->stride))
#endif

#ifndef INDEX_OF_MATRIX
#define INDEX_OF_MATRIX(a, i, j) ((i) + (j) * (a->lda))
#endif

#ifndef POINTER_OF_MATRIX
#define POINTER_OF_MATRIX(a, i, j) ((a->data) + (i) + (j) * (a->lda))
#endif

#endif /* CLINALG_MACROS_H_ */
