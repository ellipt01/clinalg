/*
 * c_linalg_cholesky.c
 *
 *  Created on: 2014/04/03
 *      Author: utsugi
 */

#include <clinalg.h>

#include "private.h"

/* cholesky decomposition */
int
c_linalg_cholesky_decomp (c_matrix *a)
{
	int			info;
	int			j;
	int			inc = 1;
	c_vector	*z;

	if (c_matrix_is_empty (a)) c_error ("c_linalg_cholesky_decomp", "matrix is empty.");
	if (!c_matrix_is_square (a)) c_error ("c_linalg_cholesky_decomp", "matrix must be square.");

	info = c_linalg_lapack_dpotrf ('U', a);

	z = c_vector_alloc (a->size1);
	c_vector_set_zero (z);
	for (j = 0; j < a->size2 - 1; j++) {
		int		n = (int) z->size - (j + 1);
		F77CALL (dcopy) (&n, z->data, &inc, POINTER_OF_MATRIX (a, j + 1, j), &inc);
	}
	c_vector_free (z);

	return info;
}

/* solve b = a * x by dpotrf and dpotrs */
int
c_linalg_cholesky_solve (c_matrix *a, c_vector *b)
{
	int			info;
	c_matrix	*c;

	if (c_matrix_is_empty (a)) c_error ("c_linalg_cholesky_decomp", "matrix is empty.");
	if (!c_matrix_is_square (a)) c_error ("c_linalg_cholesky_decomp", "matrix must be square.");

	info = c_linalg_lapack_dpotrf ('U', a);
	if (info != 0) c_error ("c_linalg_cholesky_solve", "dpotrf failed.");

	c = c_matrix_view_array (b->size, 1, b->size, b->data);
	info = c_linalg_lapack_dpotrs ('U', a, c);
	c_matrix_free (c);

	return info;
}

/* solve b = l' * l * x */
int
c_linalg_cholesky_svx (c_matrix *l, c_vector *b)
{
	int			info;
	c_matrix	*c;
	if (c_matrix_is_empty (l)) c_error ("c_linalg_cholesky_svx", "first matrix is empty.");
	if (c_vector_is_empty (b)) c_error ("c_linalg_cholesky_svx", "second vector is empty.");
	if (!c_matrix_is_square (l)) c_error ("c_linalg_cholesky_svx", "first matrix must be square.");
	if (l->size1 != b->size) c_error ("c_linalg_cholesky_svx", "matrix size does not match.");

	c = c_matrix_view_array (b->size, 1, b->size, b->data);
	info = c_linalg_lapack_dpotrs ('U', l, c);
	c_matrix_free (c);
	return info;
}

/* inverse matrix */
int
c_linalg_cholesky_invert (c_matrix *l)
{
	if (c_matrix_is_empty (l)) c_error ("c_linalg_cholesky_invert", "matrix is empty.");
	if (!c_matrix_is_square (l)) c_error ("c_linalg_cholesky_invert", "matrix must be square.");

	return c_linalg_lapack_dpotri ('U', l);
}

/* chol1up */
void
c_linalg_cholesky_1up (c_matrix *l, c_vector *u)
{
	int		n;
	int		ldr;
	double	*w;

	if (c_matrix_is_empty (l)) c_error ("cl_linalg_cholesky_1up", "matrix *r is empty");
	if (c_vector_is_empty (u)) c_error ("cl_linalg_cholesky_1up", "vector *u is empty");
	if (!c_matrix_is_square (l)) c_error ("cl_linalg_cholesky_1up", "matrix *r must be square");
	if (u->size != l->size1) c_error ("cl_linalg_cholesky_1up", "size of matrix *r and vector *u invalid");

	n = l->size1;
	ldr = l->lda;
	w = (double *) malloc (u->size * sizeof (double));
	F77CALL (dch1up) (&n, l->data, &ldr, u->data, w);
	free (w);

	return;
}

/* chol1down */
int
c_linalg_cholesky_1down (c_matrix *l, c_vector *u)
{
	int		info;
	int		n;
	int		ldr;
	double	*w;

	if (c_matrix_is_empty (l)) c_error ("cl_linalg_cholesky_1down", "matrix *r is empty");
	if (c_vector_is_empty (u)) c_error ("cl_linalg_cholesky_1down", "vector *u is empty");
	if (!c_matrix_is_square (l)) c_error ("cl_linalg_cholesky_1down", "matrix *r must be square");
	if (u->size != l->size1) c_error ("cl_linalg_cholesky_1down", "size of matrix *r and vector *u invalid");

	n = l->size1;
	ldr = l->lda;
	w = (double *) malloc (u->size * sizeof (double));
	F77CALL (dch1dn) (&n, l->data, &ldr, u->data, w, &info);
	free (w);

	return info;
}

/* cholinsert */
int
c_linalg_cholesky_insert (c_matrix *l, const int index, c_vector *u)
{
	int			n;
	int			ldr;
	int			j;
	double		*w;
	int			info;

	if (c_matrix_is_empty (l)) c_error ("c_linalg_cholesky_insert", "matrix is empty.");
	if (c_vector_is_empty (u)) c_error ("c_linalg_cholesky_insert", "vector is empty.");
	if (!c_matrix_is_square (l)) c_error ("c_linalg_cholesky_insert", "matrix must be square.");
	if (u->size != l->size1 + 1) c_error ("c_linalg_cholesky_insert", "matrix and vector size does not match.");
	if (index < 0 || l->size1 < index) c_error ("c_linalg_cholesky_insert", "index must be in [0, r->size1].");

	j = index + 1;

	n = l->size1;

	c_matrix_add_rowcols (l, 1, 1);

	ldr = l->lda;
	w = (double *) malloc (ldr * sizeof (double));
	F77CALL (dchinx) (&n, l->data, &ldr, &j, u->data, w, &info);
	free (w);

	return info;
}

/* choldelete */
void
c_linalg_cholesky_delete (c_matrix *l, const int index)
{
	int			n;
	int			ldr;
	int			j;
	double		*w;

	if (c_matrix_is_empty (l)) c_error ("c_linalg_cholesky_delete", "matrix is empty.");
	if (!c_matrix_is_square (l)) c_error ("c_linalg_cholesky_delete", "matrix must be square.");
	if (index < 0 || l->size1 <= index) c_error ("c_linalg_cholesky_delete", "index must be in [0, r->size1).");

	j = index + 1;

	n = l->size1;
	ldr = l->lda;
	w = (double *) malloc (ldr * sizeof (double));

	F77CALL (dchdex) (&n, l->data, &ldr, &j, w);
	free (w);

	c_matrix_remove_rowcols (l, 1, 1);

	return;
}
