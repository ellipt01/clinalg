/*
 * c_linalg_sv.c
 *
 *  Created on: 2014/04/19
 *      Author: utsugi
 */

#include <clinalg.h>

/* c_linalg_utils.c */
extern void	c_error (const char * function_name, const char *error_msg);

static int
_c_linalg_SV_decomp (c_matrix *a, c_matrix **u, c_matrix **vt, c_vector **s)
{
	int		info;
	char	jobu = 'A';
	char	jobvt = 'A';
	if (u == NULL) jobu = 'N';
	if (vt == NULL) jobvt = 'N';
	info = c_linalg_lapack_dgesvd (jobu, jobvt, a, u, vt, s);
	return info;
}

static int
_c_linalg_SV_decomp_econ (c_matrix *a, c_matrix **u, c_matrix **vt, c_vector **s)
{
	int		info;
	char	jobu = 'S';
	char	jobvt = 'S';
	if (u == NULL) jobu = 'N';
	if (vt == NULL) jobvt = 'N';
	info = c_linalg_lapack_dgesvd (jobu, jobvt, a, u, vt, s);
	return info;
}

int
c_linalg_SV_decomp (c_matrix *a, c_matrix **u, c_matrix **vt, c_vector **s, bool econ)
{
	int		info;
	if (c_matrix_is_empty(a)) c_error ("c_linalg_SV_decomp", "matrix *a is empty.");
	info = (econ) ? _c_linalg_SV_decomp_econ (a, u, vt, s) : _c_linalg_SV_decomp (a, u, vt, s);
	return info;
}

int
c_linalg_SV_solve (double rcond, c_matrix *a, c_vector *b, c_vector **s, int *rank)
{
	int			info;
	int			_rank;
	c_vector	*_s;

	if (c_vector_is_empty (b)) c_error ("c_linalg_SV_solve", "vector is empty.");
	if (c_matrix_is_empty (a)) c_error ("c_linalg_SV_solve", "matrix is empty.");
	if (b->stride != 1) c_error ("c_linalg_SV_solve", "cannot tread vector with stride.");
	if (a->size1 != b->size) c_error ("c_linalg_SV_solve", "vector and matrix size dose not match.");
	if (b->size < a->size2) c_vector_realloc (a->size2, b, b->size);

	{
		c_matrix	*x = c_matrix_view_array (b->size, 1, b->size, b->data);
		info = c_linalg_lapack_dgelss (rcond, a, x, &_s, &_rank);
		c_matrix_free (x);
	}

	if (info == 0 && b->size != a->size2) b->size = a->size2;

	if (rank) *rank = _rank;

	if (s) *s = _s;
	else if (!c_vector_is_empty (_s)) c_vector_free (_s);

	return info;
}

int
c_linalg_SV_lsd_solve (double rcond, c_matrix *a, c_vector *b, c_vector **s, int *rank)
{
	int			info;
	int			_rank;
	c_vector	*_s;

	if (c_vector_is_empty (b)) c_error ("c_linalg_SV_lsd_solve", "vector is empty.");
	if (c_matrix_is_empty (a)) c_error ("c_linalg_SV_lsd_solve", "matrix is empty.");
	if (b->stride != 1) c_error ("c_linalg_SV_lsd_solve", "cannot tread vector with stride.");
	if (a->size1 != b->size) c_error ("c_linalg_SV_lsd_solve", "vector and matrix size dose not match.");
	if (b->size < a->size2) c_vector_realloc (a->size2, b, b->size);

	{
		c_matrix	*x = c_matrix_view_array (b->size, 1, b->size, b->data);
		info = c_linalg_lapack_dgelsd (rcond, a, x, &_s, &_rank);
		c_matrix_free (x);
	}

	if (info == 0 && b->size != a->size2) b->size = a->size2;

	if (rank) *rank = _rank;

	if (s) *s = _s;
	else if (!c_vector_is_empty (_s)) c_vector_free (_s);

	return info;
}
