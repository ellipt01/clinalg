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

	dtrtrs_ (&uplo, &trans, &diag, &n, &nrhs, a->data, &lda, b->data, &ldb, &info);
	return info;
}

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

	dgetrf_ (&m, &n, a->data, &lda, _p->data, &info);

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
	dgesv_ (&n, &nrhs, a->data, &lda, _p->data, b->data, &ldb, &info);

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
	dgetrs_ (&trans, &n, &nrhs, lu->data, &lda, p->data, b->data, &ldb, &info);

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
	dgetri_(&n, lu->data, &lda, p->data, &wkopt, &lwork, &info);

	lwork = (int) wkopt;
	if (info != 0 || lwork <= 0) c_error ("c_linalg_lapack_dgetri", "failed to query workspace.");
	work = (double *) malloc (lwork * sizeof (double));
	if (!work) c_error ("c_linalg_lapack_dgetri", "cannot allocate memory work.");

	dgetri_(&n, lu->data, &lda, p->data, work, &lwork, &info);
	free (work);

	return info;
}

/* interface of lapack dpotrf_ */
int
c_linalg_lapack_dpotrf (char uplo, c_matrix *a)
{
	int		info;
	int		n;
	int		lda;

	if (c_matrix_is_empty (a)) c_error ("c_linalg_lapack_dpotrf", "matrix is empty.");
	if (!c_matrix_is_square (a)) c_error ("c_linalg_lapack_dpotrf", "matrix must be square.");
	if (uplo != 'L' && uplo != 'U') c_error ("c_linalg_lapack_dpotrf", "uplo must be 'L' or 'U'.");

	n = (int) a->size1;
	lda = (int) a->lda;
	dpotrf_ (&uplo, &n, a->data, &lda, &info);
	return info;
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
	return info;
}

/* interface of lapack dpotri_ */
int
c_linalg_lapack_dpotri (char uplo, c_matrix *l)
{
	int		info;
	int		n;
	int		lda;
	if (c_matrix_is_empty (l)) c_error ("c_linalg_lapack_dpotri", "matrix is empty.");
	if (!c_matrix_is_square (l)) c_error ("c_linalg_lapack_dpotri", "matrix must be square.");
	if (uplo != 'L' && uplo != 'U') c_error ("c_linalg_lapack_dpotri", "uplo must be 'L' or 'U'.");

	n = (int) l->size1;
	lda = (int)l->lda;
	dpotri_ (&uplo, &n, l->data, &lda, &info);
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
 	dgeqrf_ (&m, &n, a->data, &lda, _tau->data, &wkopt, &lwork, &info);

 	lwork = (int) wkopt;
	if (info != 0 || lwork <= 0) c_error ("c_linalg_lapack_dgeqrf", "failed to query size of workspace.");
	if ((work = (double *) malloc (lwork * sizeof (double))) == NULL)
		c_error ("c_linalg_lapack_dgeqrf", "cannot allocate memory for workspace.");
	dgeqrf_ (&m, &n, a->data, &lda, _tau->data, work, &lwork, &info);
	free (work);

	if (tau) *tau = _tau;
	else if (!c_vector_is_empty (_tau)) c_vector_free (_tau);

	return info;
}

int
c_linalg_lapack_dgeqp3 (c_matrix *a, c_vector **tau, c_vector_int **p)
{
	int			i;
	int   		info;
	int			m;
	int			n;
	int			lda;
	int			min_mn;

	double		wkopt;
	int			lwork;
	double		*work;

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
	dgeqp3_ (&m, &n, a->data, &lda, _p->data, _tau->data, &wkopt, &lwork, &info);

	lwork = (int) wkopt;
	if (info != 0 || lwork <= 0) c_error ("c_linalg_lapack_dgeqp3", "failed to query size of workspace.");
	if ((work = (double *) malloc (lwork * sizeof (double))) == NULL)
		c_error ("c_linalg_lapack_dgeqp3", "cannot allocate memory for workspace.");
	dgeqp3_ (&m, &n, a->data, &lda, _p->data, _tau->data, work, &lwork, &info);
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
	dorgqr_ (&m, &min_mn, &k, qr->data, &lda, tau->data, &wkopt, &lwork, &info);
	lwork = (int) wkopt;
	if (info != 0 || lwork <= 0) c_error ("c_linalg_lapack_dorgqr", "failed to query size of workspace.");
	if ((work = (double *) malloc (lwork * sizeof (double))) == NULL)
		c_error ("c_linalg_lapack_dorgqr", "cannot allocate memory for workspace.");
	dorgqr_ (&m, &min_mn, &k, qr->data, &lda, tau->data, work, &lwork, &info);
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
	dgels_ (&trans, &m, &n, &nrhs, a->data, &lda, b->data, &ldb, &wkopt, &lwork, &info);
	lwork = (int) wkopt;
	if (info != 0 || lwork <= 0) c_error ("c_linalg_lapack_dgels", "failed to query size of workspace.");
	if ((work = (double *) malloc (lwork * sizeof (double))) == NULL)
		c_error ("c_linalg_lapack_dgels", "cannot allocate memory for workspace.");
	dgels_ (&trans, &m, &n, &nrhs, a->data, &lda, b->data, &ldb, work, &lwork, &info);
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
	dgelsy_ (&m, &n, &nrhs, a->data, &lda, b->data, &ldb, _p->data, &rcond, &_rank, &wkopt, &lwork, &info);

	lwork = (int) wkopt;
	if (info != 0 || lwork <= 0) c_error ("c_linalg_lapack_dgelsy", "failed to query workspace");
	if ((work = (double *) malloc (lwork * sizeof (double))) == NULL)
		c_error ("c_linalg_lapack_dgelsy", "cannot allocate memory work");
	dgelsy_ (&m, &n, &nrhs, a->data, &lda, b->data, &ldb, _p->data, &rcond, &_rank, work, &lwork, &info);
	free (work);

	if (p) *p = _p;
	else if (!c_vector_int_is_empty (_p)) c_vector_int_free (_p);

	if (rank) *rank = _rank;

	return info;
}

int
c_linalg_lapack_dgesvd (char jobu, char jobvt, c_matrix *a, c_matrix **u, c_matrix **vt, c_vector **s)
{
	int			info;
	int			m;
	int			n;
	int		min_mn;
	int			lda;
	int			ldu = 1;
	double		*u_data = NULL;
	int			ldvt = 1;
	double		*vt_data = NULL;
	int			lwork;
	double		wkopt;
	double		*work;

	c_vector	*_s = NULL;
	c_matrix	*_u = NULL;
	c_matrix	*_vt = NULL;

	if (c_matrix_is_empty (a)) c_error ("c_linalg_lapack_dgesvd", "matrix is empty.");

	m = (int) a->size1;
	n = (int) a->size2;
	lda  = (int) a->lda;
	min_mn = C_MIN (a->size1, a->size2);
	_s = c_vector_alloc (min_mn);

	switch (jobu) {
		case 'A':
		case 'a':
			_u = c_matrix_alloc (a->size1, a->size1);
			u_data = _u->data;
			ldu = (int) C_MAX (1, m);
			break;

		case 'S':
		case 's':
			if (a->size1 > a->size2) _u = c_matrix_alloc (a->size1, a->size2);
			else _u = c_matrix_alloc (a->size1, a->size1);
			u_data = _u->data;
			ldu = (int) C_MAX (1, m);
			break;

		case 'O':
		case 'o':
		case 'N':
		case 'n':
			break;

		default:
			c_error ("c_linalg_lapack_dgesvd", "jobu must be 'A', 'S' 'O' or 'N'.");
			return -1;
	}

	switch (jobvt) {
		case 'A':
		case 'a':
			_vt = c_matrix_alloc (a->size2, a->size2);
			vt_data = _vt->data;
			ldvt = (int) C_MAX (1, n);
			break;

		case 'S':
		case 's':
			if (a->size1 > a->size2) _vt = c_matrix_alloc (a->size2, a->size2);
			else _vt = c_matrix_alloc (a->size1, a->size2);
			vt_data = _vt->data;
			ldvt = (int) min_mn;
			break;

		case 'O':
		case 'o':
		case 'N':
		case 'n':
			break;

		default:
			c_error ("c_linalg_lapack_dgesvd", "jobvt must be 'A', 'S', 'O' or 'N'.");
			return -1;
	}

	lwork = -1;
	dgesvd_ (&jobu, &jobvt, &m, &n, a->data, &lda, _s->data, u_data, &ldu, vt_data, &ldvt, &wkopt, &lwork, &info);

	lwork = (int) wkopt;
	if (info != 0 || lwork <= 0) c_error ("c_linalg_lapack_dgesvd", "failed to query workspace.");
	if ((work = (double *) malloc (lwork * sizeof (double))) == NULL)
		c_error ("c_linalg_lapack_dgesvd", "cannot allocate memory for work.");
	dgesvd_ (&jobu, &jobvt, &m, &n, a->data, &lda, _s->data, u_data, &ldu, vt_data, &ldvt, work, &lwork, &info);
	free (work);

	if (_s && s) *s = _s;
	if (_u && u) *u = _u;
	if (_vt && vt) *vt = _vt;

	return info;
}

int
c_linalg_lapack_dgesdd (char jobz, c_matrix *a, c_matrix **u, c_matrix **vt, c_vector **s)
{
	int			info;
	int			m;
	int			n;
	int		min_mn;
	int			lda;
	int			ldu = 1;
	double		*u_data = NULL;
	int			ldvt = 1;
	double		*vt_data = NULL;
	int			lwork;
	double		wkopt;
	double		*work;
	int			*iwork;
	c_vector	*_s = NULL;
	c_matrix	*_u = NULL;
	c_matrix	*_vt = NULL;

	if (c_matrix_is_empty (a)) c_error ("c_linalg_lapack_dgesdd", "matrix *a is empty.");

	m = (int) a->size1;
	n = (int) a->size2;
	lda  = (int) a->lda;
	min_mn = C_MIN (a->size1, a->size2);
	iwork = (int *) malloc (8 * min_mn * sizeof (int));
	_s = c_vector_alloc (min_mn);

	switch (jobz) {
		case 'A':
		case 'a':
			_u = c_matrix_alloc (a->size1, a->size1);
			_vt = c_matrix_alloc (a->size2, a->size2);
			u_data  = _u->data;
			vt_data = _vt->data;
			ldu  = (int) C_MAX (1, m);
			ldvt = (int) C_MAX (1, n);
			break;

		case 'S':
		case 's':
			if (a->size1 > a->size2) {
				_u = c_matrix_alloc (a->size1, a->size2);
				_vt = c_matrix_alloc (a->size2, a->size2);
			} else {
				_u = c_matrix_alloc (a->size1, a->size1);
				_vt = c_matrix_alloc (a->size1, a->size2);
			}
			ldu = (int) C_MAX (1, m);
			ldvt = (int) min_mn;
			u_data  = _u->data;
			vt_data = _vt->data;
			break;

		case 'N':
		case 'n':
			break;

		default:
			c_error ("c_linalg_lapack_dgesdd", "jobz must be 'A' 'S' or 'N'.");
			return -1;
	}

	lwork = -1;
	dgesdd_ (&jobz, &m, &n, a->data, &lda, _s->data, u_data, &ldu, vt_data, &ldvt, &wkopt, &lwork, iwork, &info);

	lwork = (int) wkopt;
	if (info != 0 || lwork <= 0) c_error ("c_linalg_lapack_dgesdd", "failed to query workspace.");
	if ((work = (double *) malloc (lwork * sizeof (double))) == NULL)
		c_error ("c_linalg_lapack_dgesdd", "failed to allocate memory for work.");

	dgesdd_ (&jobz, &m, &n, a->data, &lda, _s->data, u_data, &ldu, vt_data, &ldvt, work, &lwork, iwork, &info);
	free (work);
	free (iwork);

	if (_s && s) *s = _s;
	if (_u && u) *u = _u;
	if (_vt && vt) *vt = _vt;

	return info;
}

int
c_linalg_lapack_dgelss (double rcond, c_matrix *a, c_matrix *b, c_vector **s, int *rank)
{
	int			info;
	int			m;
	int			n;
	int		min_mn;
	int			nrhs;
	int			lda;
	int			ldb;
	int			lwork;
	double		wkopt;
	double		*work;
	int			_rank;
	c_vector	*_s;

	if (c_matrix_is_empty(a)) c_error ("c_linalg_lapack_dgelss", "matrix *a is empty.");
	if (c_matrix_is_empty(b)) c_error ("c_linalg_lapack_dgelss", "matrix *b is empty.");

	m = (int) a->size1;
	n = (int) a->size2;
	nrhs = (int) b->size2;
	lda = (int) a->lda;
	ldb = (int) C_MAX (1, (int) C_MAX (m, n));
	min_mn = C_MIN (m, n);

	_s = c_vector_alloc (min_mn);

	lwork = -1;
	dgelss_ (&m, &n, &nrhs, a->data, &lda, b->data, &ldb, _s->data, &rcond, &_rank, &wkopt, &lwork, &info);
	lwork = (int) wkopt;
	if (info != 0 || lwork <= 0) c_error ("c_linalg_lapack_dgelss", "failed to query workspace.");
	if ((work = (double *) malloc (lwork * sizeof (double))) == NULL)
		c_error ("c_linalg_lapack_dgelss", "failed to allocate memory.");
	dgelss_ (&m, &n, &nrhs, a->data, &lda, b->data, &ldb, _s->data, &rcond, &_rank, work, &lwork, &info);
	free (work);

	if (rank) *rank = (int) _rank;

	if (s) *s = _s;
	else if (!c_vector_is_empty (_s)) c_vector_free (_s);

	return info;
}

int
c_linalg_lapack_dgelsd (double rcond, c_matrix *a, c_matrix *b, c_vector **s, int *rank)
{
	int			info;
	int			m;
	int			n;
	int			min_mn;
	int			nrhs;
	int			lda;
	int			ldb;
	int			lwork;
	double		wkopt;
	double		*work;
	int			smlsiz;
	int			nlvl;
	int			liwork;
	int			*iwork;
	int			_rank;
	c_vector	*_s;

	if (c_matrix_is_empty(a)) c_error ("c_linalg_lapack_dgelsd", "matrix *a is empty.");
	if (c_matrix_is_empty(b)) c_error ("c_linalg_lapack_dgelsd", "matrix *b is empty.");

	m = (int) a->size1;
	n = (int) a->size2;
	nrhs = (int) b->size2;
	lda = (int) a->lda;
	ldb = (int) C_MAX (1, C_MAX (m, n));
	min_mn = C_MIN (m, n);

	/* iwork */
	{
		int	ispec = 9;
		int	n1 = 0;
		int	n2 = 0;
		int	n3 = 0;
		int	n4 = 0;
		smlsiz = ilaenv_ (&ispec, "DGELSD", " ", &n1, &n2, &n3, &n4);
	}
	nlvl = (int) C_MAX (0, (int) (log2 ((double) min_mn / (double) (smlsiz + 1))) + 1);
	liwork = (int) C_MAX (1, 3 * min_mn * nlvl + 11 * min_mn);
	if ((iwork = (int *) malloc (liwork * sizeof (int))) == NULL)
		c_error ("c_linalg_lapack_dgelsd", "failed to allocate iwork.");

	_s = c_vector_alloc (min_mn);

	lwork = -1;
	dgelsd_ (&m, &n, &nrhs, a->data, &lda, b->data, &ldb, _s->data, &rcond, &_rank, &wkopt, &lwork, iwork, &info);
	lwork = (int) wkopt;

	if (info != 0 || lwork <= 0) c_error ("c_linalg_lapack_dgelsd", "failed to query workspace.");
	if ((work = (double *) malloc (((int) lwork) * sizeof (double))) == NULL)
		c_error ("c_linalg_lapack_dgelsd", "allocation of work failed.");
	dgelsd_ (&m, &n, &nrhs, a->data, &lda, b->data, &ldb, _s->data, &rcond, &_rank, work, &lwork, iwork, &info);
	free (work);
	free (iwork);

	if (rank) *rank = (int) _rank;

	if (s) *s = _s;
	else if (!c_vector_is_empty (_s)) c_vector_free (_s);

	return info;
}
