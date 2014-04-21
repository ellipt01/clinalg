/*
 * c_linalg_sv.c
 *
 *  Created on: 2014/04/19
 *      Author: utsugi
 */

#include <math.h>
#include <clinalg.h>

/* c_linalg_util.c */
extern void	c_error (const char * function_name, const char *error_msg);

/* lapack */
#ifndef HAVE_LAPACK_H
extern int		ilaenv_ (int *ispec, char *name, char *opts, int *n1, int *n2, int *n3, int *n4);
extern void	dgesvd_ (char *jobu, char *jobvt, int *m, int *n, double *a_data, int *lda, double *s_data,
						  double *u_data, int *ldu, double *vt_data, int *ldvt, double *w, int *lwork, int *info);
extern void	dgesdd_ (char *jobz, int *m, int *n, double *a_data, int *lda, double *s_data, double *u_data, int *ldu,
						  double *vt_data, int *ldvt, double *work, int *lwork, int *iwork, int *info);
extern void	dgelss_ (int *m, int *n, int *nrhs, double *a_data, int *lda, double *b_data, int *ldb,
						  double *s_data, double *rcond, int *lrank, double *w, int *lwork, int *info);
extern void	dgelsd_ (int *m, int *n, int *nrhs, double *a_data, int *lda, double *b_data, int *ldb,
						  double *s_data, double *rcond, int *lrank, double *w, int *lwork, int *iwork, int *info);
#endif

int
c_linalg_lapack_dgesvd (char jobu, char jobvt, c_matrix *a, c_matrix **u, c_matrix **vt, c_vector **s)
{
	int			info;
	int			m;
	int			n;
	size_t		min_mn;
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
	size_t		min_mn;
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
	size_t		min_mn;
	int			nrhs;
	int			lda;
	int			ldb;
	int			lwork;
	double		wkopt;
	double		*work;
	int			_rank;
	c_vector	*_s;

	if (c_matrix_is_empty(a)) c_error ("c_linalg_lapack_dgelss", "input matrix *a is empty.");
	if (c_matrix_is_empty(b)) c_error ("c_linalg_lapack_dgelss", "input matrix *b is empty.");

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
	size_t		min_mn;
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

	if (c_matrix_is_empty(a)) c_error ("c_linalg_lapack_dgelsd", "input matrix *a is empty.");
	if (c_matrix_is_empty(b)) c_error ("c_linalg_lapack_dgelsd", "input matrix *b is empty.");

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

	_s = c_vector_alloc ((size_t) min_mn);

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

static int
_c_linalg_SV_decomp (c_matrix *a, c_matrix **u, c_matrix **vt, c_vector **s)
{
	int		info;
	char	jobu = 'A';
	char	jobvt = 'A';
	if (u == NULL) jobu = 'N';
	if (vt == NULL) jobvt = 'N';
	info = c_linalg_lapack_dgesvd (jobu, jobvt, a, u, vt, s);
	return info;
}

static int
_c_linalg_SV_decomp_econ (c_matrix *a, c_matrix **u, c_matrix **vt, c_vector **s)
{
	int		info;
	char	jobu = 'S';
	char	jobvt = 'S';
	if (u == NULL) jobu = 'N';
	if (vt == NULL) jobvt = 'N';
	info = c_linalg_lapack_dgesvd (jobu, jobvt, a, u, vt, s);
	return info;
}

int
c_linalg_SV_decomp (c_matrix *a, c_matrix **u, c_matrix **vt, c_vector **s, bool econ)
{
	int		info;
	if (c_matrix_is_empty(a)) c_error ("c_linalg_SV_decomp", "matrix *a is empty.");
	info = (econ) ? _c_linalg_SV_decomp_econ (a, u, vt, s) : _c_linalg_SV_decomp (a, u, vt, s);
	return info;
}

int
c_linalg_SV_solve (double rcond, c_matrix *a, c_vector *b, c_vector **s, int *rank)
{
	int			info;
	int			_rank;
	c_vector	*_s;

	if (c_vector_is_empty (b)) c_error ("c_linalg_SV_solve", "vector is empty.");
	if (c_matrix_is_empty (a)) c_error ("c_linalg_SV_solve", "matrix is empty.");
	if (a->size1 != b->size) c_error ("c_linalg_SV_solve", "vector and matrix size dose not match.");
	if (b->size < a->size2) b->data = (double *) realloc (b->data, a->size2 * sizeof (double));

	{
		c_matrix	*x = c_matrix_view_array (a->size1, 1, a->lda, b->data);
		info = c_linalg_lapack_dgelss (rcond, a, x, &_s, &_rank);
		c_matrix_free (x);
	}

	if (b->size != a->size2) b->size = a->size2;

	if (rank) *rank = _rank;

	if (s) *s = _s;
	else if (!c_vector_is_empty (_s)) c_vector_free (_s);

	return info;
}

int
c_linalg_SV_lsd_solve (double rcond, c_matrix *a, c_vector *b, c_vector **s, int *rank)
{
	int			info;
	int			_rank;
	c_vector	*_s;

	if (c_vector_is_empty (b)) c_error ("c_linalg_SV_lsd_solve", "vector is empty.");
	if (c_matrix_is_empty (a)) c_error ("c_linalg_SV_lsd_solve", "matrix is empty.");
	if (a->size1 != b->size) c_error ("c_linalg_SV_lsd_solve", "vector and matrix size dose not match.");
	if (b->size != a->size2) b->data = (double *) realloc (b->data, a->size2 * sizeof (double));

	{
		c_matrix	*x = c_matrix_view_array (a->size1, 1, a->lda, b->data);
		info = c_linalg_lapack_dgelsd (rcond, a, x, &_s, &_rank);
		c_matrix_free (x);
	}

	if (b->size != a->size2) b->size = a->size2;

	if (rank) *rank = _rank;

	if (s) *s = _s;
	else if (!c_vector_is_empty (_s)) c_vector_free (_s);

	return info;
}
