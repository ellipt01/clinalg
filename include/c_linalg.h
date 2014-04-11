/*
 * c_linalg.h
 *
 *  Wrapper of lapack, blas and qrupdate with a GSL like interface
 *
 *  Created on: 2014/04/08
 *      Author: utsugi
 */

#ifndef C_LINALG_H_
#define C_LINALG_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <c_matrix.h>
#include <c_linalg_cholesky.h>

#define C_MAX(a, b)		((a >= b) ? a : b)
#define C_MIN(a, b)		((a <= b) ? a : b)

#define DBL_EPSILON		2.2204460492503131e-16

#ifdef __cplusplus
}
#endif

#endif /* C_LINALG_H_ */
