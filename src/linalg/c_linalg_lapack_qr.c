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
c_linalg_lapack_dtrtrs (char uplo, char trans, char diag, c_matrix *a, c_matrix *b)
{
	int		n;
	int		nrhs;
	int		lda;
	int		ldb;
	int		info;

	if (uplo != 'U' && uplo != 'u' && uplo != 'l' && uplo != 'L') c_error ("c_linalg_lapack_dtrtrs", "uplo must be 'U' or 'L'.");
	if (trans != 'T' && trans != 't' && trans != 'N' && trans != 'n') c_error ("c_linalg_lapack_dtrtrs", "trans must be 'T' or 'N'.");
	if (diag != 'U' && diag != 'u' && diag != 'N' && diag != 'n') c_error ("c_linalg_lapack_dtrtrs", "diag must be 'U' or 'N'.");
	if (c_matrix_is_empty (a)) c_error ("c_linalg_lapack_dtrtrs", "matrix *a is empty.");
	if (c_matrix_is_empty (b)) c_error ("c_linalg_lapack_dtrtrs", "matrix *b is empty.");

	n = (int) a->size1;
	nrhs = (int) b->size2;
	lda = (int) a->lda;
	ldb = (int) b->lda;

	F77CALL (dtrtrs) (&uplo, &trans, &diag, &n, &nrhs, a->data, &lda, b->data, &ldb, &info);
	return info;
}

int
c_linalg_lapack_dgeqrf (c_matrix *a, c_vector **tau)
{
	int   		info;
	int			m;
	int			n;
	int			lda;

	int			min_mn;
	int			ltau;

	double		wkopt;
	int			lwork;
	double		*work;

	c_vector	*_tau;

	if (c_matrix_is_empty (a)) c_error ("c_linalg_lapack_dgeqrf", "matrix is empty.");

	m = (int) a->size1;
	n = (int) a->size2;
	lda = (int) a->lda;

	min_mn = C_MIN (a->size1, a->size2);
	ltau = C_MAX (min_mn, 1);

	_tau = c_vector_alloc (ltau);

 	lwork = -1;
 	F77CALL (dgeqrf) (&m, &n, a->data, &lda, _tau->data, &wkopt, &lwork, &info);

 	lwork = (int) wkopt;
	if (info != 0 || lwork <= 0) c_error ("c_linalg_lapack_dgeqrf", "failed to query size of workspace.");
	if ((work = (double *) malloc (lwork * sizeof (double))) == NULL)
		c_error ("c_linalg_lapack_dgeqrf", "cannot allocate memory for workspace.");
	F77CALL (dgeqrf) (&m, &n, a->data, &lda, _tau->data, work, &lwork, &info);
	free (work);

	if (tau) *tau = _tau;
	else if (!c_vector_is_empty (_tau)) c_vector_free (_tau);

	return info;
}

int
c_linalg_lapack_dgeqp3 (c_matrix *a, c_vector **tau, c_vector_int **p)
{
	int				i;
	int   			info;
	int				m;
	int				n;
	int				lda;
	int				min_mn;

	double			wkopt;
	int				lwork;
	double			*work;

	c_vector		*_tau;
	c_vector_int	*_p;

	if (c_matrix_is_empty (a)) c_error ("c_linalg_lapack_dgeqp3", "matrix is empty.");

	m = (int) a->size1;
	n = (int) a->size2;
	lda = (int) a->lda;
	min_mn = (int) C_MIN (a->size1, a->size2);

	_tau = c_vector_alloc (min_mn);
	_p = c_vector_int_alloc (n);
	for (i = 0; i < _p->size; i++) _p->data[i] = i + 1;

	lwork = -1;
	F77CALL (dgeqp3) (&m, &n, a->data, &lda, _p->data, _tau->data, &wkopt, &lwork, &info);

	lwork = (int) wkopt;
	if (info != 0 || lwork <= 0) c_error ("c_linalg_lapack_dgeqp3", "failed to query size of workspace.");
	if ((work = (double *) malloc (lwork * sizeof (double))) == NULL)
		c_error ("c_linalg_lapack_dgeqp3", "cannot allocate memory for workspace.");
	F77CALL (dgeqp3) (&m, &n, a->data, &lda, _p->data, _tau->data, work, &lwork, &info);
	free (work);

	if (p) *p = _p;
	else if (!c_vector_int_is_empty (_p)) c_vector_int_free (_p);

	if (tau) *tau = _tau;
	else if (!c_vector_is_empty (_tau)) c_vector_free (_tau);

	return info;
}

int
c_linalg_lapack_dorgqr (c_matrix *qr, const c_vector *tau)
{
	int 		info;
	int			m;
	int			n;
	int			min_mn;
	int			k;
	int			lda;
	int			lwork;
	double		wkopt;
	double		*work;

	if (c_matrix_is_empty (qr)) c_error ("c_linalg_lapack_dorgqr", "matrix is empty.");
	if (c_vector_is_empty (tau)) c_error ("c_linalg_lapack_dorgqr", "vector *tau is empty.");

	m = (int) qr->size1;
	n = (int) qr->size2;
	min_mn = (int) C_MIN (m, n);
	k = (int) tau->size;
	lda = (int) qr->lda;

 	lwork = -1;
 	F77CALL (dorgqr) (&m, &min_mn, &k, qr->data, &lda, tau->data, &wkopt, &lwork, &info);
	lwork = (int) wkopt;
	if (info != 0 || lwork <= 0) c_error ("c_linalg_lapack_dorgqr", "failed to query size of workspace.");
	if ((work = (double *) malloc (lwork * sizeof (double))) == NULL)
		c_error ("c_linalg_lapack_dorgqr", "cannot allocate memory for workspace.");
	F77CALL (dorgqr) (&m, &min_mn, &k, qr->data, &lda, tau->data, work, &lwork, &info);
	free (work);

	return info;
}

int
c_linalg_lapack_dgels (char trans, c_matrix *a, c_matrix *b)
{
	int			info;
	int			m;
	int			n;
	int			nrhs;
	int			lda;
	int			ldb;
	double		wkopt;
	double		*work;
	int			lwork;

	if (c_matrix_is_empty(a)) c_error ("c_linalg_lapack_dgels", "matrix *a is empty.");
	if (c_matrix_is_empty(b)) c_error ("c_linalg_lapack_dgels", "matrix *b is empty.");

	m = (int) a->size1;
	n = (int) a->size2;
	nrhs = (int) b->size2;
	lda = (int) a->lda;
	ldb = (int) C_MAX (1, C_MAX (m, n));

	lwork = -1;
	F77CALL (dgels) (&trans, &m, &n, &nrhs, a->data, &lda, b->data, &ldb, &wkopt, &lwork, &info);
	lwork = (int) wkopt;
	if (info != 0 || lwork <= 0) c_error ("c_linalg_lapack_dgels", "failed to query size of workspace.");
	if ((work = (double *) malloc (lwork * sizeof (double))) == NULL)
		c_error ("c_linalg_lapack_dgels", "cannot allocate memory for workspace.");
	F77CALL (dgels) (&trans, &m, &n, &nrhs, a->data, &lda, b->data, &ldb, work, &lwork, &info);
	free (work);

	return info;
}

int
c_linalg_lapack_dgelsy (double rcond, c_matrix *a, c_matrix *b, c_vector_int **p, int *rank)
{
	int			info;
	int			m;
	int			n;
	int			nrhs;
	int			lda;
	int			ldb;
	int			_rank;
	double		wkopt;
	double		*work;
	int			lwork;

	c_vector_int	*_p;

	if (c_matrix_is_empty (a)) c_error ("c_linalg_lapack_dgelsy", "matrix *a is empty");
	if (c_matrix_is_empty (b)) c_error ("c_linalg_lapack_dgelsy", "matrix *b is empty");

	m = (int) a->size1;
	n = (int) a->size2;
	nrhs = (int) b->size2;
	lda = (int) a->lda;
	ldb = (int) C_MAX (1, C_MAX (m, n));
	_p = c_vector_int_alloc (a->size2);

	lwork = -1;
	F77CALL (dgelsy) (&m, &n, &nrhs, a->data, &lda, b->data, &ldb, _p->data, &rcond, &_rank, &wkopt, &lwork, &info);

	lwork = (int) wkopt;
	if (info != 0 || lwork <= 0) c_error ("c_linalg_lapack_dgelsy", "failed to query workspace");
	if ((work = (double *) malloc (lwork * sizeof (double))) == NULL)
		c_error ("c_linalg_lapack_dgelsy", "cannot allocate memory work");
	F77CALL (dgelsy) (&m, &n, &nrhs, a->data, &lda, b->data, &ldb, _p->data, &rcond, &_rank, work, &lwork, &info);
	free (work);

	if (p) *p = _p;
	else if (!c_vector_int_is_empty (_p)) c_vector_int_free (_p);

	if (rank) *rank = _rank;

	return info;
}

