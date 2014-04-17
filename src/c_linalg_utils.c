/*
 * c_linalg_utils.c
 *
 *  Created on: 2014/04/10
 *      Author: utsugi
 */

#include <c_linalg.h>

void
c_error (const char * function_name, const char *error_msg)
{
	fprintf (stderr, "ERROR: %s: %s\n", function_name, error_msg);
	exit (1);
}

c_matrix *
c_linalg_permutation_matrix_row (const size_t size1, const size_t size2, const c_vector_int *p)
{
	int			i;
	size_t		min;
	c_matrix	*d;

	if (c_vector_int_is_empty (p)) c_error ("c_linalg_permutation_matrix_row", "vector_int is empty.");

	d = c_matrix_identity (size1, size2);
	min = C_MIN (size1, p->size);
	for (i = 0; i < min; i++) {
		int		pi = c_vector_int_get (p, i) - 1;
		if (pi < size1 && pi != i) c_matrix_swap_rows (i, pi, d);
	}
	return d;
}

c_matrix *
c_linalg_permutation_matrix_col (const size_t size1, const size_t size2, const c_vector_int *p)
{
	int			i;
	size_t		min;
	c_matrix	*d;

	if (c_vector_int_is_empty (p)) c_error ("c_linalg_permutation_matrix_row", "vector_int is empty.");

	d = c_matrix_identity (size1, size2);
	min = C_MIN (size2, p->size);
	for (i = 0; i < min; i++) {
		int		pi = c_vector_int_get (p, i) - 1;
		if (pi < size1 && pi != i) c_matrix_swap_cols (i, pi, d);
	}
	return d;
}
