/*
 * c_linalg_lapack_eigen.c
 *
 *  Created on: 2015/09/01
 *      Author: utsugi
 */

#include <clinalg.h>

#include "private.h"

/* real symmetric matrix */
int
c_linalg_lapack_dsyev (const char jobz, const char uplo, c_matrix *a, c_vector **w)
{
	int			info;
	int			n;
	int			lda;
	double		wkopt;
	double		*work;
	int			lwork;

	c_vector	*_w;

	if (c_matrix_is_empty (a)) c_error ("c_linalg_lapack_dsyev", "matrix A is empty");
	if (!c_matrix_is_square (a)) c_error ("c_linalg_lapack_dsyev", "matrix A must be square");
	if (jobz != 'N' && jobz != 'V') c_error ("c_linalg_lapack_dsyev", "jobz must be 'N' or 'V'");
	if (uplo != 'U' && uplo != 'L') c_error ("c_linalg_lapack_dsyev", "uplo must be 'U' or 'L'");

	n = (int) a->size1;
	lda = (int) a->lda;

	_w = c_vector_alloc (n);

	lwork = -1;
	F77CALL (dsyev) (&jobz, &uplo, &n, a->data, &lda, _w->data, &wkopt, &lwork, &info);

	lwork = (int) wkopt;
	if ((int) info != 0 || lwork <= 0) c_error ("c_linalg_lapack_dsyev", "failed to query workspace");
	if ((work = (double *) malloc (lwork * sizeof (double))) == NULL)
		c_error ("c_linalg_lapack_dsyev", "failed to allocate memory");
	F77CALL (dsyev) (&jobz, &uplo, &n, a->data, &lda, _w->data, work, &lwork, &info);

	if (w) *w = _w;
	else c_vector_free (_w);

	free (work);

	return info;
}

/* real symmetric matrix */
int
c_linalg_lapack_dsyevd (const char jobz, const char uplo, c_matrix *a, c_vector **w)
{
	int			info;
	int			n;
	int			lda;
	double		wkopt;
	double		*work;
	int			lwork;
	int			*iwork;
	int			liwork;

	c_vector	*_w;

	if (c_matrix_is_empty (a)) c_error ("c_linalg_lapack_dsyevd", "matrix A is empty");
	if (!c_matrix_is_square (a)) c_error ("c_linalg_lapack_dsyevd", "matrix A must be square");
	if (jobz != 'N' && jobz != 'V') c_error ("c_linalg_lapack_dsyevd", "jobz must be 'N' or 'V'");
	if (uplo != 'U' && uplo != 'L') c_error ("c_linalg_lapack_dsyevd", "uplo must be 'U' or 'L'");

	n = (int) a->size1;
	lda = (int) a->lda;

	_w = c_vector_alloc (n);

	liwork = (jobz == 'N') ? 1 : 3 + 5 * n;
	if ((iwork = (int *) malloc (liwork * sizeof (int))) == NULL)
		c_error ("c_linalg_lapack_dsyevd", "failed to allocate memory");

	lwork  = -1;
	F77CALL (dsyevd) (&jobz, &uplo, &n, a->data, &lda, _w->data, &wkopt, &lwork, iwork, &liwork, &info);

	lwork = (int) wkopt;
	if ((int) info != 0 || lwork <= 0) c_error ("c_linalg_lapack_dsyevd", "failed to query workspace");
	if ((work = (double *) malloc (lwork * sizeof (double))) == NULL)
		c_error ("c_linalg_lapack_dsyevd", "failed to allocate memory");

	F77CALL (dsyevd) (&jobz, &uplo, &n, a->data, &lda, _w->data, work, &lwork, iwork, &liwork, &info);

	if (w) *w = _w;
	else c_vector_free (_w);

	free (work);
	free (iwork);

	return info;
}

/* real unsymmetric matrix */
int
c_linalg_lapack_dgeev (const char jobvl, const char jobvr, c_matrix *a, c_vector **wr, c_vector **wi, c_matrix **vl, c_matrix **vr)
{
	int			info;

	int			n;
	int			lda;
	double		*__wr;
	double		*__wi;
	double		*__vl;
	int			ldvl;
	double		*__vr;
	int			ldvr;
	double		wkopt;
	double		*work;
	int			lwork;

	c_vector	*_wr;
	c_vector	*_wi;

	c_matrix	*_vl = NULL;
	c_matrix	*_vr = NULL;

	if (c_matrix_is_empty (a)) c_error ("c_linalg_lapack_dgeev", "matrix A is empty");

	n = (int) a->size1;
	lda = (int) a->lda;

	switch (jobvl) {
		case 'V':
		case 'v':
			ldvl = n;
			_vl = c_matrix_alloc ((size_t) ldvl, (size_t) n);
			__vl = _vl->data;
			break;

		case 'N':
		case 'n':
			ldvl = 1;
			__vl = NULL;
			break;

		default:
			c_error ("c_linalg_lapack_dgeev", "jobvl must be 'V' or 'N'");
			return -1;
	}

	switch (jobvr) {
		case 'V':
		case 'v':
			ldvr = n;
			_vr = c_matrix_alloc ((size_t) ldvr, (size_t) n);
			__vr = _vr->data;
			break;

		case 'N':
		case 'n':
			ldvr = 1;
			__vr = NULL;
			break;

		default:
			c_error ("c_linalg_lapack_dgeev", "jobvr must be 'V' or 'N'");
			return -1;
	}

	_wr = c_vector_alloc (n);
	__wr = _wr->data;
	_wi = c_vector_alloc (n);
	__wi = _wi->data;

	lwork = -1;
	F77CALL (dgeev) (&jobvl, &jobvr, &n, a->data, &lda, __wr, __wi, __vl, &ldvl, __vr, &ldvr, &wkopt, &lwork, &info);

	lwork = (int) wkopt;
	if ((int) info != 0 || lwork <= 0) c_error ("c_linalg_lapack_dgeev", "failed to query workspace");
	if ((work = (double *) malloc (lwork * sizeof (double))) == NULL)
		c_error ("c_linalg_lapack_dgeev", "failed to allocate memory");
	F77CALL (dgeev) (&jobvl, &jobvr, &n, a->data, &lda, __wr, __wi, __vl, &ldvl, __vr, &ldvr, work, &lwork, &info);

	if (vl) *vl = _vl;
	else if (_vl) c_matrix_free (_vl);

	if (vr && _vr) *vr = _vr;
	else if (_vr) c_matrix_free (_vr);

	if (wr) *wr = _wr;
	else c_vector_free (_wr);

	if (wi) *wi = _wi;
	else c_vector_free (_wi);

	free (work);

	return info;
}

/* schur decomposition */
int
c_linalg_lapack_dgees (const char jobvs, const char sort, c_matrix *a, c_vector **wr, c_vector **wi, c_matrix **vs, void *select)
{
	int			n;
	int			lda;
	int			sdim;
	double		*work;
	int			*bwork;
	double		wkopt;
	int			lwork;

	int			ldvs = 1;
	c_vector	*_wr = NULL;
	double		*__wr = NULL;
	c_vector	*_wi = NULL;
	double		*__wi = NULL;
	c_matrix	*_vs = NULL;
	double		*__vs = NULL;

	int			info;

	if (c_matrix_is_empty (a)) c_error ("c_linalg_clapack_dgees", "matrix a is empty");
	if (!c_matrix_is_square (a)) c_error ("c_linalg_clapack_dgees", "matrix a must be square");

	n = (int) a->size1;
	lda = (int) a->lda;

	switch (sort) {
		case 'S':
		case 's':
			if (select == NULL) c_error ("c_linalg_clapack_dgees", "void *select is NULL");
		break;

		case 'N':
		case 'n':
		default:
		break;
	}

	switch (jobvs) {
		case 'V':
		case 'v':
			ldvs = n;
			_vs = c_matrix_alloc (a->size1, a->size2);
			__vs = _vs->data;
		break;

		case 'N':
		case 'n':
		default:
		break;
	}

	_wr = c_vector_alloc (n);
	__wr = _wr->data;
	_wi = c_vector_alloc (n);
	__wi = _wi->data;
	bwork = (int *) malloc (n * sizeof (int));

	lwork = -1;
	F77CALL (dgees) (&jobvs, &sort, select, &n, a->data, &lda, &sdim, __wr, __wi, __vs, &ldvs, &wkopt, &lwork, bwork, &info);
	lwork = (int) wkopt;
	work = (double *) malloc (lwork * sizeof (double));
	F77CALL (dgees) (&jobvs, &sort, select, &n, a->data, &lda, &sdim, __wr, __wi, __vs, &ldvs, work, &lwork, bwork, &info);

	if (wr) *wr = _wr;
	else c_vector_free (_wr);

	if (wi) *wi = _wi;
	else c_vector_free (_wi);

	if (vs) *vs = _vs;
	else c_matrix_free (_vs);

	free (work);
	free (bwork);

	return info;
}
