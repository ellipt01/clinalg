/*
 * test_matrix.c
 *
 *  Created on: 2014/10/01
 *      Author: utsugi
 */

#include <math.h>
#include <clinalg.h>

#include "test_clinalg.h"

extern int		size1;
extern int		size2;

c_matrix *
create_symmetric_matrix (int n)
{
	c_matrix	*tmp = random_matrix (n, n);
	c_matrix	*c = c_matrix_dot_matrix_transpose (1., tmp, tmp);
	c_matrix_free (tmp);

	return c;
}

bool
test_matrix_mv (void)
{
	c_matrix	*a = create_symmetric_matrix (size1);
	c_vector	*x = random_vector (size1);

	c_vector	*v0 = c_matrix_dot_vector (1., a, x);
	c_vector	*v1 = c_matrix_symm_upper_dot_vector (1., a, x);
	c_vector	*v2 = c_matrix_symm_lower_dot_vector (1., a, x);

	/* v1 = -v0 + v1 */
	c_vector_axpy (-1., v0, v1);
	/* v2 = -v0 + v2 */
	c_vector_axpy (-1., v0, v2);
	return (c_vector_nrm (v1) < 1.e-3 && c_vector_nrm (v2) < 1.e-3);
}

bool
test_matrix_mm (void)
{
	c_matrix	*a = create_symmetric_matrix (size1);
	c_matrix	*b = random_matrix (size1, size2);
	c_matrix	*m0 = c_matrix_dot_matrix (1., a, b);
	c_matrix	*m1 = c_matrix_symm_upper_dot_matrix (1., a, b);
	c_matrix	*m2 = c_matrix_symm_lower_dot_matrix (1., a, b);

	c_matrix_sub (m1, m0);
	c_matrix_sub (m2, m0);
	return (c_matrix_nrm (m1, '1') < 1.e-3 && c_matrix_nrm (m2, '1') < 1.e-3);
}
