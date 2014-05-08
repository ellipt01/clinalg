/*
 * test_clinalg.h
 *
 *  Created on: 2014/04/14
 *      Author: utsugi
 */

#ifndef TEST_CLINALG_H_
#define TEST_CLINALG_H_

c_vector	*random_vector (size_t size);
c_matrix	*random_matrix (size_t size1, size_t size2);

bool	test_cholesky_decomp (void);
bool	test_cholesky_svx (void);
bool	test_cholesky_1up (void);
bool	test_cholesky_1down (void);
bool	test_cholesky_insert (void);
bool	test_cholesky_delete (void);

bool	test_LU_decomp (void);
bool	test_LU_solve (void);
bool	test_LU_svx (void);
bool	test_LU_invert (void);
bool	test_LU_1up (void);

bool	test_QR_decomp (void);
bool	test_QR_decomp_econ (void);
bool	test_QR_solve (void);
bool	test_lsQ_solve (void);
bool	test_QR_Rsolve (void);
bool	test_QR_1up (void);
bool	test_QR_1down (void);
bool	test_QR_colinsert (void);
bool	test_QR_coldelete (void);
bool	test_QR_rowinsert (void);
bool	test_QR_rowdelete (void);

bool	test_SV_decomp (void);
bool	test_SV_solve (void);
bool	test_SV_lsd_solve (void);

#endif /* TEST_CLINALG_H_ */
