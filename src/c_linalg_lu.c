/*
 * c_linalg_lu.c
 *
 *  Created on: 2014/04/03
 *      Author: utsugi
 */

#include <c_matrix.h>

/* c_linalg_util.c */
extern void	c_error (const char * function_name, const char *error_msg);

/* lapack */
extern void	dgetrf_ (long *m, long *n, double *data, long *lda, long *ipiv, long *info);
extern void	dgetri_ (long *n, double *data, long *lda, long *ipiv, double *work, long *lwork, long *info);
extern void	dgetrs_ (char *trans, long *n, long *nrhs, double *a, long *lda, long *ipiv, double *b, long *ldb, long *info);

int
c_linalg_lapack_dgetrf (c_matrix *a, long **p)
{
	long		info;
 	long		m;
 	long		n;
 	long		lda;
 	long		*_p;

	if (c_matrix_is_empty (a)) c_error ("c_linalg_lapack_dgetrf", "matrix is empty.");

 	m = (long) a->size1;
 	n = (long) a->size2;
 	lda  = (long) a->lda;
 	_p = (long *) malloc (a->size1 * sizeof (long));
	dgetrf_ (&m, &n, a->data, &lda, _p, &info);
	if (p) *p = _p;

	return (int) info;
}

int
c_linalg_lapack_dgetrs (char trans, c_matrix *lu, c_matrix *b, long *p)
{
	long	info;
	long	n;
	long	nrhs;
	long	lda;
	long	ldb;

	if (c_matrix_is_empty (lu)) c_error ("c_linalg_lapack_dgetrs", "matrix LU is empty");
	if (!c_matrix_is_square (lu)) c_error ("c_linalg_lapack_dgetrs", "matrix LU must be square");
	if (!p) c_error ("c_linalg_lapack_dgetrs", "permutation is empty");
	if (trans != 'N' && trans != 'T' && trans != 'C') c_error ("c_linalg_lapack_dgetrs", "trans must be 'N', 'T' or 'C'");

	n = (long) lu->size1;
	nrhs = (long) b->size2;
	lda = (long) lu->lda;
	ldb = (long) b->lda;
	dgetrs_ (&trans, &n, &nrhs, lu->data, &lda, p, b->data, &ldb, &info);

	return (int) info;
}

int
c_linalg_lapack_dgetri (c_matrix *a, long *p)
{
	long		info;
 	long		n;
	long		lda;
	double		*work;
	double		wkopt;
	long		lwork;

	if (c_matrix_is_empty (a)) c_error ("c_linalg_lapack_dgetri", "matrix *a is empty.");
	if (!c_matrix_is_square (a)) c_error ("c_linalg_lapack_dgetri", "matrix *a must be square.");
	if (!p) c_error ("c_linalg_lapack_dgetri", "permutation is empty");

	n = (long) a->size1;
	lda  = (long) a->lda;

	lwork = -1;
	dgetri_(&n, a->data, &lda, p, &wkopt, &lwork, &info);

	lwork = (long) wkopt;
	if ((int) info != 0 || lwork <= 0) c_error ("c_linalg_lapack_dgetri", "failed to query workspace.");
	work = (double *) malloc (lwork * sizeof (double));
	if (!work) c_error ("c_linalg_lapack_dgetri", "cannot allocate memory work.");

	dgetri_(&n, a->data, &lda, p, work, &lwork, &info);
	free (work);

	return (int) info;
}

int
c_linalg_LU_decomp (c_matrix *a, long **p)
{
	long		info;
	long		*_p;

	if (c_matrix_is_empty (a)) c_error ("c_linalg_LU_decomp", "matrix is empty.");

	info = c_linalg_lapack_dgetrf (a, &_p);
	if (p) *p = _p;

	return info;
}

int
c_linalg_LU_solve (c_matrix *a, c_vector *b, long *p)
{
	int  		info;
	c_matrix	*c;

	if (c_matrix_is_empty (a)) c_error ("c_linalg_LU_solve", "matrix is empty.");
	if (c_vector_is_empty (b)) c_error ("c_linalg_LU_solve", "vector is empty.");
	if (!p) c_error ("c_linalg_LU_solve", "permutation is empty.");
	if (b->size != a->size1) c_error ("c_linalg_LU_solve", "matrix and vector size dose not match.");
	if (!c_matrix_is_square (a)) c_error ("c_linalg_LU_solve", "matrix must be square.");

	c = c_matrix_view_array (b->size, 1, b->size, b->data);
	info = c_linalg_lapack_dgetrs ('N', a, c, p);
	c_matrix_free (c);

	return info;
}

int
c_linalg_LU_invert (c_matrix *a, long *p)
{
	int			info;

	if (c_matrix_is_empty (a)) c_error ("c_linalg_LU_invert", "matrix is empty.");
	if (!c_matrix_is_square (a)) c_error ("c_linalg_LU_invert", "matrix must be square.");

	info = c_linalg_lapack_dgetri (a, p);

	return (int) info;
}
