/*
 * c_linalg.c
 *
 *  Created on: 2014/05/13
 *      Author: utsugi
 */

#include <math.h>
#include <clinalg.h>

#include "private.h"

int
c_linalg_lapack_dgetrf (c_matrix *a, c_vector_int **p)
{
	int				i;
	int				info;
 	int				m;
 	int				n;
 	int				lda;
 	c_vector_int	*_p;

 	if (c_matrix_is_empty (a)) c_error ("c_linalg_lapack_dgetrf", "matrix is empty.");

 	m = (int) a->size1;
 	n = (int) a->size2;
 	lda  = (int) a->lda;
 	_p = c_vector_int_alloc (a->size1);
	for (i = 0; i < _p->size; i++) _p->data[i] = i + 1;

	F77CALL (dgetrf) (&m, &n, a->data, &lda, _p->data, &info);

	if (p) *p = _p;
	else if (!c_vector_int_is_empty (_p)) c_vector_int_free (_p);

	return info;
}

int
c_linalg_lapack_dgesv (c_matrix *a, c_matrix *b, c_vector_int **p)
{
	int				info;
	int				n;
	int				nrhs;
	int				lda;
	int				ldb;
	c_vector_int	*_p;

	if (c_matrix_is_empty (a)) c_error ("c_linalg_lapack_dgesv", "matrix is empty.");
	if (!c_matrix_is_square (a)) c_error ("c_linalg_lapack_dgesv", "matrix must be square.");

	n = (int) a->size1;
	nrhs = (int) b->size2;
	lda = (int) a->lda;
	ldb = (int) b->lda;

	_p = c_vector_int_alloc (n);
	F77CALL (dgesv) (&n, &nrhs, a->data, &lda, _p->data, b->data, &ldb, &info);

	if (p) *p = _p;
	else if (!c_vector_int_is_empty (_p)) c_vector_int_free (_p);

	return info;
}

int
c_linalg_lapack_dgetrs (char trans, c_matrix *lu, c_matrix *b, c_vector_int *p)
{
	int			info;
	int			n;
	int			nrhs;
	int			lda;
	int			ldb;

	if (c_matrix_is_empty (lu)) c_error ("c_linalg_lapack_dgetrs", "matrix is empty.");
	if (!c_matrix_is_square (lu)) c_error ("c_linalg_lapack_dgetrs", "matrix must be square.");
	if (!p) c_error ("c_linalg_lapack_dgetrs", "permutation is empty.");
	if (trans != 'N' && trans != 'T' && trans != 'C') c_error ("c_linalg_lapack_dgetrs", "trans must be 'N', 'T' or 'C'.");

	n = (int) lu->size1;
	nrhs = (int) b->size2;
	lda = (int) lu->lda;
	ldb = (int) b->lda;
	F77CALL (dgetrs) (&trans, &n, &nrhs, lu->data, &lda, p->data, b->data, &ldb, &info);

	return info;
}

int
c_linalg_lapack_dgetri (c_matrix *lu, c_vector_int *p)
{
	int			info;
 	int			n;
	int			lda;
	double		*work;
	double		wkopt;
	int			lwork;

	if (c_matrix_is_empty (lu)) c_error ("c_linalg_lapack_dgetri", "matrix is empty.");
	if (!c_matrix_is_square (lu)) c_error ("c_linalg_lapack_dgetri", "matrix must be square.");
	if (!p) c_error ("c_linalg_lapack_dgetri", "permutation is empty.");

	n = (int) lu->size1;
	lda  = (int) lu->lda;

	lwork = -1;
	F77CALL (dgetri) (&n, lu->data, &lda, p->data, &wkopt, &lwork, &info);

	lwork = (int) wkopt;
	if (info != 0 || lwork <= 0) c_error ("c_linalg_lapack_dgetri", "failed to query workspace.");
	work = (double *) malloc (lwork * sizeof (double));
	if (!work) c_error ("c_linalg_lapack_dgetri", "cannot allocate memory work.");

	F77CALL (dgetri) (&n, lu->data, &lda, p->data, work, &lwork, &info);
	free (work);

	return info;
}

