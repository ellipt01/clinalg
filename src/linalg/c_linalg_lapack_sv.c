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
c_linalg_lapack_dgesvd (char jobu, char jobvt, c_matrix *a, c_matrix **u, c_matrix **vt, c_vector **s)
{
	int			info;
	int			m;
	int			n;
	int			min_mn;
	int			lda;
	int			ldu;
	double		*__u;
	int			ldvt;
	double		*__vt;
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
			__u = _u->data;
			ldu = (int) C_MAX (1, m);
			break;

		case 'S':
		case 's':
			if (a->size1 > a->size2) _u = c_matrix_alloc (a->size1, a->size2);
			else _u = c_matrix_alloc (a->size1, a->size1);
			__u = _u->data;
			ldu = (int) C_MAX (1, m);
			break;

		case 'O':
		case 'o':
		case 'N':
		case 'n':
			ldu = 1;
			__u = NULL;
			break;

		default:
			c_error ("c_linalg_lapack_dgesvd", "jobu must be 'A', 'S' 'O' or 'N'.");
			return -1;
	}

	switch (jobvt) {
		case 'A':
		case 'a':
			_vt = c_matrix_alloc (a->size2, a->size2);
			__vt = _vt->data;
			ldvt = (int) C_MAX (1, n);
			break;

		case 'S':
		case 's':
			if (a->size1 > a->size2) _vt = c_matrix_alloc (a->size2, a->size2);
			else _vt = c_matrix_alloc (a->size1, a->size2);
			__vt = _vt->data;
			ldvt = (int) min_mn;
			break;

		case 'O':
		case 'o':
		case 'N':
		case 'n':
			ldvt = 1;
			__vt = NULL;
			break;

		default:
			c_error ("c_linalg_lapack_dgesvd", "jobvt must be 'A', 'S', 'O' or 'N'.");
			return -1;
	}

	lwork = -1;
	F77CALL (dgesvd) (&jobu, &jobvt, &m, &n, a->data, &lda, _s->data, __u, &ldu, __vt, &ldvt, &wkopt, &lwork, &info);

	lwork = (int) wkopt;
	if (info != 0 || lwork <= 0) c_error ("c_linalg_lapack_dgesvd", "failed to query workspace.");
	if ((work = (double *) malloc (lwork * sizeof (double))) == NULL)
		c_error ("c_linalg_lapack_dgesvd", "cannot allocate memory for work.");
	F77CALL (dgesvd) (&jobu, &jobvt, &m, &n, a->data, &lda, _s->data, __u, &ldu, __vt, &ldvt, work, &lwork, &info);
	free (work);

	if (s) *s = _s;
	else if (!c_vector_is_empty (_s)) c_vector_free (_s);

	if (u) *u = _u;
	else if (_u) c_matrix_free (_u);

	if (vt) *vt = _vt;
	else if (_vt) c_matrix_free (_vt);

	return info;
}

int
c_linalg_lapack_dgesdd (char jobz, c_matrix *a, c_matrix **u, c_matrix **vt, c_vector **s)
{
	int			info;
	int			m;
	int			n;
	int			min_mn;
	int			lda;
	int			ldu;
	double		*__u;
	int			ldvt;
	double		*__vt;
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
			__u  = _u->data;
			__vt = _vt->data;
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
			__u  = _u->data;
			__vt = _vt->data;
			break;

		case 'N':
		case 'n':
			ldu = 1;
			ldvt = 1;
			__u = NULL;
			__vt = NULL;
			break;

		default:
			c_error ("c_linalg_lapack_dgesdd", "jobz must be 'A' 'S' or 'N'.");
			return -1;
	}

	lwork = -1;
	F77CALL (dgesdd) (&jobz, &m, &n, a->data, &lda, _s->data, __u, &ldu, __vt, &ldvt, &wkopt, &lwork, iwork, &info);

	lwork = (int) wkopt;
	if (info != 0 || lwork <= 0) c_error ("c_linalg_lapack_dgesdd", "failed to query workspace.");
	if ((work = (double *) malloc (lwork * sizeof (double))) == NULL)
		c_error ("c_linalg_lapack_dgesdd", "failed to allocate memory for work.");

	F77CALL (dgesdd) (&jobz, &m, &n, a->data, &lda, _s->data, __u, &ldu, __vt, &ldvt, work, &lwork, iwork, &info);
	free (work);
	free (iwork);

	if (s) *s = _s;
	else if (!c_vector_is_empty (_s)) c_vector_free (_s);

	if (u) *u = _u;
	else if (_u) c_matrix_free (_u);

	if (vt) *vt = _vt;
	else if (_vt) c_matrix_free (_vt);

	return info;
}

int
c_linalg_lapack_dgelss (double rcond, c_matrix *a, c_matrix *b, c_vector **s, int *rank)
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
	F77CALL (dgelss) (&m, &n, &nrhs, a->data, &lda, b->data, &ldb, _s->data, &rcond, &_rank, &wkopt, &lwork, &info);
	lwork = (int) wkopt;
	if (info != 0 || lwork <= 0) c_error ("c_linalg_lapack_dgelss", "failed to query workspace.");
	if ((work = (double *) malloc (lwork * sizeof (double))) == NULL)
		c_error ("c_linalg_lapack_dgelss", "failed to allocate memory.");
	F77CALL (dgelss) (&m, &n, &nrhs, a->data, &lda, b->data, &ldb, _s->data, &rcond, &_rank, work, &lwork, &info);
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
		smlsiz = F77CALL (ilaenv) (&ispec, "DGELSD", " ", &n1, &n2, &n3, &n4);
	}
	nlvl = (int) C_MAX (0, (int) (log2 ((double) min_mn / (double) (smlsiz + 1))) + 1);
	liwork = (int) C_MAX (1, 3 * min_mn * nlvl + 11 * min_mn);
	if ((iwork = (int *) malloc (liwork * sizeof (int))) == NULL)
		c_error ("c_linalg_lapack_dgelsd", "failed to allocate iwork.");

	_s = c_vector_alloc (min_mn);

	lwork = -1;
	F77CALL (dgelsd) (&m, &n, &nrhs, a->data, &lda, b->data, &ldb, _s->data, &rcond, &_rank, &wkopt, &lwork, iwork, &info);
	lwork = (int) wkopt;

	if (info != 0 || lwork <= 0) c_error ("c_linalg_lapack_dgelsd", "failed to query workspace.");
	if ((work = (double *) malloc (((int) lwork) * sizeof (double))) == NULL)
		c_error ("c_linalg_lapack_dgelsd", "allocation of work failed.");
	F77CALL (dgelsd) (&m, &n, &nrhs, a->data, &lda, b->data, &ldb, _s->data, &rcond, &_rank, work, &lwork, iwork, &info);
	free (work);
	free (iwork);

	if (rank) *rank = (int) _rank;

	if (s) *s = _s;
	else if (!c_vector_is_empty (_s)) c_vector_free (_s);

	return info;
}
