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
extern void	dpotrf_ (char *uplo, long *n, double *a, long *lda, long *info);
extern void	dpotrs_ (char *uplo, long *n, long *nrhs, double *a, long *lda, double *b, long *ldb, long *info);
extern void	dpotri_ (char *uplo, long *n, double *a, long *lda, long *info);
/* qrupdate: cholinsert/delete */
extern void	dchinx_ (int *n, double *R, int *ldr, int *j, double *u, double *w, int *info);
extern void	dchdex_ (int *n, double *R, int *ldr, int *j, double *w);

/* interface of lapack dpotrf_ */
int
c_linalg_lapack_dpotrf (char uplo, c_matrix *a)
{
	long	info;
	long	n;
	long	lda;

	if (c_matrix_is_empty (a)) c_error ("c_linalg_lapack_dpotrf", "matrix is empty.");
	if (!c_matrix_is_square (a)) c_error ("c_linalg_lapack_dpotrf", "matrix must be square.");
	if (uplo != 'L' && uplo != 'U') c_error ("c_linalg_lapack_dpotrf", "uplo must be 'L' or 'U'.");

	n = (long) a->size1;
	lda = (long) a->lda;
	dpotrf_ (&uplo, &n, a->data, &lda, &info);
	return (int) info;
}

/* interface of lapack dpotrs_ */
int
c_linalg_lapack_dpotrs (char uplo, c_matrix *a, c_matrix *b)
{
	long		info;
	long		n;
	long		nrhs;
	long		lda;
	long		ldb;

	if (c_matrix_is_empty (a)) c_error ("c_linalg_lapack_dpotrs", "first matrix is empty.");
	if (c_matrix_is_empty (b)) c_error ("c_linalg_lapack_dpotrs", "second matrix is empty.");
	if (!c_matrix_is_square (a)) c_error ("c_linalg_lapack_dpotrs", "first matrix must be square.");
	if (a->size1 != b->size1) c_error ("c_linalg_lapack_dpotrs", "matrix size does not match.");
	if (uplo != 'L' && uplo != 'U') c_error ("c_linalg_lapack_dpotrs", "uplo must be 'L' or 'U'");

	n = (long) a->size1;
	nrhs = (long) b->size2;
	lda = (long) a->lda;
	ldb = (long) b->lda;
	dpotrs_ (&uplo, &n, &nrhs, a->data, &lda, b->data, &ldb, &info);
	return (int) info;
}

/* interface of lapack dpotri_ */
int
c_linalg_lapack_dpotri (char uplo, c_matrix *a)
{
	long	info;
	long	n;
	long	lda;
	if (c_matrix_is_empty (a)) c_error ("c_linalg_lapack_dpotri", "matirx is empty");
	if (!c_matrix_is_square (a)) c_error ("c_linalg_lapack_dpotri", "matirx must be square");
	if (uplo != 'L' && uplo != 'U') c_error ("c_linalg_lapack_dpotri", "uplo must be 'L' or 'U'");

	n = (long) a->size1;
	lda = (long) a->lda;
	dpotri_ (&uplo, &n, a->data, &lda, &info);
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
c_linalg_cholesky_svx (c_matrix *a, c_vector *b)
{
	int			info;
	c_matrix	*c;
	if (c_matrix_is_empty (a)) c_error ("c_linalg_cholesky_svx", "first matrix is empty.");
	if (c_vector_is_empty (b)) c_error ("c_linalg_cholesky_svx", "second vector is empty.");
	if (!c_matrix_is_square (a)) c_error ("c_linalg_cholesky_svx", "first matrix must be square.");
	if (a->size1 != b->size) c_error ("c_linalg_cholesky_svx", "matrix size does not match.");

	c = c_matrix_view_array (b->size, 1, b->size, b->data);
	info = c_linalg_lapack_dpotrs ('U', a, c);
	c_matrix_free (c);
	return info;
}

/* inverse matrix */
int
c_linalg_cholesky_invert (c_matrix *a)
{
	int			info;

	if (c_matrix_is_empty (a)) c_error ("c_linalg_cholesky_invert", "matrix is empty.");
	if (!c_matrix_is_square (a)) c_error ("c_linalg_cholesky_invert", "matrix must be square.");

	info = c_linalg_lapack_dpotrf ('U', a);
	if (info != 0) c_error ("c_linalg_cholesky_invert", "DPOTRF failed.");

	return c_linalg_lapack_dpotri ('U', a);
}

/* cholinsert */
int
c_linalg_cholesky_insert (c_matrix *r, const int index, const c_vector *u)
{
	int		n;
	int		ldr;
	int		j;
	double	*w;
	int		info;

	if (c_matrix_is_empty (r)) c_error ("c_linalg_cholesky_insert", "matrix is empty.");
	if (c_vector_is_empty (u)) c_error ("c_linalg_cholesky_insert", "vector is empty.");
	if (r->size1 != r->size2) c_error ("c_linalg_cholesky_insert", "matrix must be square.");
	if (u->size != r->size1 + 1) c_error ("c_linalg_cholesky_insert", "matrix and vector size does not match.");
	if (index < 0 || r->size1 < index) c_error ("c_linalg_cholesky_insert", "index must be in [0, r->size1].");

	j = index + 1;

	n = r->size1;

	c_matrix_add_col (r);
	c_matrix_add_row (r);

	ldr = r->lda;
	w = (double *) malloc (ldr * sizeof (double));
	dchinx_ (&n, r->data, &ldr, &j, u->data, w, &info);
	free (w);

	return info;
}

/* choldelete */
void
c_linalg_cholesky_delete (c_matrix *r, const int index)
{
	int		n;
	int		ldr;
	int		j;
	double	*w;

	if (c_matrix_is_empty (r)) c_error ("c_linalg_cholesky_delete", "matrix is empty.");
	if (r->size1 != r->size2) c_error ("c_linalg_cholesky_delete", "matrix must be square.");
	if (index < 0 || r->size1 <= index) c_error ("c_linalg_cholesky_delete", "index must be in [0, r->size1).");

	j = index + 1;

	n = r->size1;
	ldr = r->lda;
	w = (double *) malloc (ldr * sizeof (double));

	dchdex_ (&n, r->data, &ldr, &j, w);
	free (w);
	c_matrix_remove_col (r);
	if (!c_matrix_is_empty (r)) c_matrix_remove_row (r);

	return;
}