/*
 * c_linalg_cholesky.c
 *
 *  Created on: 2014/04/03
 *      Author: utsugi
 */

#include <c_matrix.h>

/* c_linalg_util.c */
extern void	c_error (const char * function_name, const char *error_msg);

/* lapack: cholesky decomposition */
extern void	dpotrf_ (char *uplo, int *n, double *a, int *lda, int *info);
extern void	dpotrs_ (char *uplo, int *n, int *nrhs, double *a, int *lda, double *b, int *ldb, int *info);
extern void	dpotri_ (char *uplo, int *n, double *a, int *lda, int *info);
/* qrupdate: cholinsert/delete */
extern void	dchinx_ (int *n, double *L, int *ldr, int *j, double *u, double *w, int *info);
extern void	dchdex_ (int *n, double *L, int *ldr, int *j, double *w);

/* interface of lapack dpotrf_ */
int
c_linalg_lapack_dpotrf (char uplo, c_matrix *a)
{
	int	info;
	int	n;
	int	lda;

	if (c_matrix_is_empty (a)) c_error ("c_linalg_lapack_dpotrf", "matrix is empty.");
	if (!c_matrix_is_square (a)) c_error ("c_linalg_lapack_dpotrf", "matrix must be square.");
	if (uplo != 'L' && uplo != 'U') c_error ("c_linalg_lapack_dpotrf", "uplo must be 'L' or 'U'.");

	n = (int) a->size1;
	lda = (int) a->lda;
	dpotrf_ (&uplo, &n, a->data, &lda, &info);
	return (int) info;
}

/* interface of lapack dpotrs_ */
int
c_linalg_lapack_dpotrs (char uplo, c_matrix *l, c_matrix *b)
{
	int		info;
	int		n;
	int		nrhs;
	int		lda;
	int		ldb;

	if (c_matrix_is_empty (l)) c_error ("c_linalg_lapack_dpotrs", "first matrix is empty.");
	if (c_matrix_is_empty (b)) c_error ("c_linalg_lapack_dpotrs", "second matrix is empty.");
	if (!c_matrix_is_square (l)) c_error ("c_linalg_lapack_dpotrs", "first matrix must be square.");
	if (l->size1 != b->size1) c_error ("c_linalg_lapack_dpotrs", "matrix size does not match.");
	if (uplo != 'L' && uplo != 'U') c_error ("c_linalg_lapack_dpotrs", "uplo must be 'L' or 'U'.");

	n = (int) l->size1;
	nrhs = (int) b->size2;
	lda = (int) l->lda;
	ldb = (int) b->lda;
	dpotrs_ (&uplo, &n, &nrhs, l->data, &lda, b->data, &ldb, &info);
	return (int) info;
}

/* interface of lapack dpotri_ */
int
c_linalg_lapack_dpotri (char uplo, c_matrix *l)
{
	int	info;
	int	n;
	int	lda;
	if (c_matrix_is_empty (l)) c_error ("c_linalg_lapack_dpotri", "matrix is empty.");
	if (!c_matrix_is_square (l)) c_error ("c_linalg_lapack_dpotri", "matrix must be square.");
	if (uplo != 'L' && uplo != 'U') c_error ("c_linalg_lapack_dpotri", "uplo must be 'L' or 'U'.");

	n = (int) l->size1;
	lda = (int)l->lda;
	dpotri_ (&uplo, &n, l->data, &lda, &info);
	return (int) info;
}

/* cholesky decomposition */
int
c_linalg_cholesky_decomp (c_matrix *a)
{
	int			info;
	if (c_matrix_is_empty (a)) c_error ("c_linalg_cholesky_decomp", "matrix is empty.");
	if (!c_matrix_is_square (a)) c_error ("c_linalg_cholesky_decomp", "matrix must be square.");

	info = c_linalg_lapack_dpotrf ('U', a);
	return info;
}

/* solve b = a * x */
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
	int			info;

	if (c_matrix_is_empty (l)) c_error ("c_linalg_cholesky_invert", "matrix is empty.");
	if (!c_matrix_is_square (l)) c_error ("c_linalg_cholesky_invert", "matrix must be square.");

	info = c_linalg_lapack_dpotrf ('U', l);
	if (info != 0) c_error ("c_linalg_cholesky_invert", "DPOTRF failed.");

	return c_linalg_lapack_dpotri ('U', l);
}

/* cholinsert */
int
c_linalg_cholesky_insert (c_matrix *l, const int index, const c_vector *u)
{
	int		n;
	int		ldr;
	int		j;
	double	*w;
	int		info;

	if (c_matrix_is_empty (l)) c_error ("c_linalg_cholesky_insert", "matrix is empty.");
	if (c_vector_is_empty (u)) c_error ("c_linalg_cholesky_insert", "vector is empty.");
	if (!c_matrix_is_square (l)) c_error ("c_linalg_cholesky_insert", "matrix must be square.");
	if (u->size != l->size1 + 1) c_error ("c_linalg_cholesky_insert", "matrix and vector size does not match.");
	if (index < 0 || l->size1 < index) c_error ("c_linalg_cholesky_insert", "index must be in [0, r->size1].");

	j = index + 1;

	n = l->size1;

	c_matrix_add_col (l);
	c_matrix_add_row (l);

	ldr = l->lda;
	w = (double *) malloc (ldr * sizeof (double));
	dchinx_ (&n, l->data, &ldr, &j, u->data, w, &info);
	free (w);

	return info;
}

/* choldelete */
void
c_linalg_cholesky_delete (c_matrix *l, const int index)
{
	int		n;
	int		ldr;
	int		j;
	double	*w;

	if (c_matrix_is_empty (l)) c_error ("c_linalg_cholesky_delete", "matrix is empty.");
	if (!c_matrix_is_square (l)) c_error ("c_linalg_cholesky_delete", "matrix must be square.");
	if (index < 0 || l->size1 <= index) c_error ("c_linalg_cholesky_delete", "index must be in [0, r->size1).");

	j = index + 1;

	n = l->size1;
	ldr = l->lda;
	w = (double *) malloc (ldr * sizeof (double));

	dchdex_ (&n, l->data, &ldr, &j, w);
	free (w);
	c_matrix_remove_col (l);
	if (!c_matrix_is_empty (l)) c_matrix_remove_row (l);

	return;
}
