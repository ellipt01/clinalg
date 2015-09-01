/*
 * c_linalg_lu.c
 *
 *  Created on: 2014/04/03
 *      Author: utsugi
 */

#include <clinalg.h>

#include "private.h"

int
c_linalg_LU_decomp (c_matrix *a, c_vector_int **p)
{
	int				info;
	c_vector_int	*_p = NULL;

	if (c_matrix_is_empty (a)) c_error ("c_linalg_LU_decomp", "matrix is empty.");

	info = c_linalg_lapack_dgetrf (a, &_p);

	if (p) *p = _p;
	else if (!c_vector_int_is_empty (_p)) c_vector_int_free (_p);

	return info;
}

void
c_linalg_LU_unpack (const c_matrix *lu, c_matrix **l, c_matrix **u)
{
	int		min_mn;

	if (c_matrix_is_empty (lu)) c_error ("c_linalg_LU_unpack", "matrix is empty.");

	min_mn = (int) C_MIN (lu->size1, lu->size2);

	if (l) {
		int			i;
		c_matrix	*_l = c_matrix_alloc (lu->size1, min_mn);
		c_matrix_set_zero (_l);
		c_matrix_lower_triangular_memcpy (_l, lu);
		for (i = 0; i < min_mn; i++) c_matrix_set (_l, i, i, 1.);
		*l = _l;
	}

	if (u) {
		c_matrix	*_u = c_matrix_alloc (min_mn, lu->size2);
		c_matrix_set_zero (_u);
		c_matrix_upper_triangular_memcpy (_u, lu);
		*u = _u;
	}
	return;
}

int
c_linalg_LU_solve (c_matrix *a, c_vector *b, c_vector_int **p)
{
	int  			info;
	c_matrix		*c;
	c_vector_int	*_p = NULL;

	if (c_matrix_is_empty (a)) c_error ("c_linalg_LU_solve", "matrix is empty.");
	if (c_vector_is_empty (b)) c_error ("c_linalg_LU_solve", "vector is empty.");
	if (!c_matrix_is_square (a)) c_error ("c_linalg_LU_solve", "matrix must be square.");
	if (b->stride != 1) c_error ("c_linalg_LU_solve", "cannot tread vector with stride.");
	if (b->size != a->size1) c_error ("c_linalg_LU_solve", "matrix and vector size dose not match.");

	c = c_matrix_view_array (b->size, 1, b->size, b->data);
	info = c_linalg_lapack_dgesv (a, c, &_p);
	c_matrix_free (c);

	if (p) *p = _p;
	else if (!c_vector_int_is_empty (_p)) c_vector_int_free (_p);

	return info;
}

int
c_linalg_LU_svx (c_matrix *lu, c_vector *b, c_vector_int *p)
{
	int  			info;
	c_matrix		*c;

	if (c_matrix_is_empty (lu)) c_error ("c_linalg_LU_svx", "matrix is empty.");
	if (c_vector_is_empty (b)) c_error ("c_linalg_LU_svx", "vector is empty.");
	if (c_vector_int_is_empty (p)) c_error ("c_linalg_LU_svx", "permutation is empty.");
	if (!c_matrix_is_square (lu)) c_error ("c_linalg_LU_svx", "matrix must be square.");
	if (b->stride != 1) c_error ("c_linalg_LU_svx", "cannot tread vector with stride.");
	if (b->size != lu->size1) c_error ("c_linalg_LU_svx", "matrix and vector size dose not match.");

	c = c_matrix_view_array (b->size, 1, b->size, b->data);
	info = c_linalg_lapack_dgetrs ('N', lu, c, p);
	c_matrix_free (c);

	return info;
}

int
c_linalg_LU_invert (c_matrix *lu, c_vector_int *p)
{
	int		info;

	if (c_matrix_is_empty (lu)) c_error ("c_linalg_LU_invert", "matrix is empty.");
	if (c_vector_int_is_empty (p)) c_error ("c_linalg_LU_invert", "permutation is empty.");
	if (!c_matrix_is_square (lu)) c_error ("c_linalg_LU_invert", "matrix must be square.");

	info = c_linalg_lapack_dgetri (lu, p);

	return info;
}

void
c_linalg_LU_1up (c_matrix *l, c_matrix *u, c_vector_int *p, c_vector *s, c_vector *t)
{
	int			m;
	int			n;
	int			ldl;
	int			ldu;
	double		*w;

	if (c_matrix_is_empty (l)) c_error ("c_linalg_LU_1up", "matrix is empty.");
	if (c_matrix_is_empty (u)) c_error ("c_linalg_LU_1up", "matrix is empty.");
	if (c_vector_is_empty (s)) c_error ("c_linalg_LU_1up", "vector *s is empty.");
	if (c_vector_is_empty (t)) c_error ("c_linalg_LU_1up", "vector *t is empty.");
	if (c_vector_int_is_empty (p)) c_error ("c_linalg_LU_1up", "permulation is empty.");
	if (s->size != l->size1) c_error ("c_linalg_LU_1up", "vector and matrix size dose not match.");
	if (t->size != u->size2) c_error ("c_linalg_LU_1up", "vector and matrix size dose not match.");
	if (s->stride != 1 || t->stride != 1) c_error ("c_linalg_LU_1up", "cannot tread vector with stride.");

	m = (int) l->size1;
	n = (int) u->size2;
	ldl = (int) l->lda;
	ldu = (int) u->lda;
	w = (double *) malloc (l->size1 * sizeof (double));
	F77CALL (dlup1up) (&m, &n, l->data, &ldl, u->data, &ldu, p->data, s->data, t->data, w);
	free (w);

	return;
}
