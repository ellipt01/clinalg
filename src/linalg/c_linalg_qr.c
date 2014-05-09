/*
 * c_linalg_qr.c
 *
 *  Created on: 2014/04/11
 *      Author: utsugi
 */

#include <clinalg.h>

/* c_linalg_util.c */
extern void	c_error (const char * function_name, const char *error_msg);

/* lapack */
#ifndef HAVE_LAPACK_H
extern void	dgeqrf_ (int *m, int *n, double *data, int *lda, double *tau, double *work, int *lwork, int *info);
extern void	dgeqp3_ (int *m, int *n, double *a, int *lda, int *jpvt, double *tau, double *work, int *lwork, int *info);
extern void	dorgqr_ (int *m, int *n, int *k, double *data, int *lda, double *tau, double *work, int *lwork, int *info);
extern void	dgels_  (char *trans, int *m, int *n, int *nrhs, double *a_data, int *lda, double *b_data, int *ldb, double *w, int *lwork, int *info);
extern void	dgelsy_ (int *m, int *n, int *nrhs, double *a, int *lda, double *b, int *ldb, int *jpvt, double *rcond, int *rank, double *work, int *lwork, int *info);
extern void	dtrsv_ (char *uplo, char *trans, char *diag, int *n, double *r, int *lda, double *y, int *incy);
#endif

/* qrupdate*/
extern void	dqr1up_ (int *m, int *n, int *k, double *Q, int *ldq, double *R, int *ldr, double *u, double *v, double *w);
extern void	dqrinc_ (int *m, int *n, int *k, double *Q, int *ldq, double *R, int *ldr, int *j, double *x, double *w);
extern void	dqrinr_ (int *m, int *n, double *Q, int *ldq, double *R, int *ldr, int *j, double *x, double *w);
extern void	dqrdec_ (int *m, int *n, int *k, double *Q, int *ldq, double *R, int *ldr, int *j, double *w);
extern void	dqrder_ (int *m, int *n, double *Q, int *ldq, double *R, int *ldr, int *j, double *w);

int
c_linalg_lapack_dgeqrf (c_matrix *a, c_vector **tau)
{
	int   		info;
	int			m;
	int			n;
	int			lda;

	size_t		min_mn;
	size_t		ltau;

	double		wkopt;
	int			lwork;
	double		*work;

	c_vector	*_tau;

	if (c_matrix_is_empty (a)) c_error ("c_linalg_lapack_dgeqrf", "matrix is empty.");

	m = (int) a->size1;
	n = (int) a->size2;
	lda = (int) a->lda;

	min_mn = (size_t) C_MIN (a->size1, a->size2);
	ltau = (size_t) C_MAX (min_mn, 1);

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
	size_t		min_mn;

	double		wkopt;
	int			lwork;
	double		*work;

	c_vector		*_tau;
	c_vector_int	*_p;

	if (c_matrix_is_empty (a)) c_error ("c_linalg_lapack_dgeqp3", "matrix is empty.");

	m = (int) a->size1;
	n = (int) a->size2;
	lda = (int) a->lda;
	min_mn = (size_t) C_MIN (a->size1, a->size2);

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
c_linalg_QR_decomp (c_matrix *a, c_vector_int **p, c_vector **tau)
{
	int				info;
 	c_vector_int	*_p = NULL;
	c_vector		*_tau = NULL;

	if (c_matrix_is_empty (a)) c_error ("c_linalg_QR_decomp", "matrix is empty.");

	if (p == NULL) info = c_linalg_lapack_dgeqrf (a, &_tau);
	else info = c_linalg_lapack_dgeqp3 (a, &_tau, &_p);

	if (p) *p = _p;
	else if (!c_vector_int_is_empty (_p)) c_vector_int_free (_p);

	if (tau) *tau = _tau;
	else if (!c_vector_is_empty (_tau)) c_vector_free (_tau);

	return info;
}

static int
_c_linalg_QR_unpack (const c_matrix *qr, const c_vector *tau, c_matrix **q, c_matrix **r)
{
	int			info = 0;
	c_matrix	*_q;
	c_matrix	*_r;

	if (c_matrix_is_empty (qr)) c_error ("c_linalg_QR_unpack", "matrix is empty.");
	if (c_vector_is_empty (tau)) c_error ("c_linalg_QR_unpack", "vector is empty.");

	if (r) {
		_r = c_matrix_alloc (qr->size1, qr->size2);
		c_matrix_set_zero (_r);
		c_matrix_upper_triangular_memcpy (_r, qr);
		*r = _r;
	}

	if (q) {
		size_t	min_mn = (size_t) C_MIN (qr->size1, qr->size2);
		_q = c_matrix_alloc (qr->size1, qr->size1);
		if (qr->size1 == qr->size2) c_matrix_memcpy (_q, qr);
		else c_matrix_mncopy (_q, 0, 0, qr->size1, min_mn, qr);
		info = c_linalg_lapack_dorgqr (_q, tau);
		*q = _q;
	}

	return info;
}

static int
_c_linalg_QR_unpack_econ (const c_matrix *qr, const c_vector *tau, c_matrix **q, c_matrix **r)
{
	int			info = 0;
	int			min_mn;
	c_matrix	*_q;
	c_matrix	*_r;

	if (c_matrix_is_empty (qr)) c_error ("c_linalg_QR_unpack", "matrix is empty.");
	if (c_vector_is_empty (tau)) c_error ("c_linalg_QR_unpack", "vector is empty.");

	min_mn = (int) C_MIN (qr->size1, qr->size2);

	if (r) {
		c_matrix	*_qr = c_matrix_submatrix (0, 0, min_mn, qr->size2, qr);
		_r = c_matrix_alloc (min_mn, qr->size2);
		c_matrix_set_zero (_r);
		c_matrix_upper_triangular_memcpy (_r, _qr);
		c_matrix_free (_qr);
		*r = _r;
	}

	if (q) {
		c_matrix	*_qr = c_matrix_submatrix (0, 0, qr->size1, min_mn, qr);
		_q = c_matrix_alloc (qr->size1, min_mn);
		c_matrix_memcpy (_q, _qr);
		c_matrix_free (_qr);
		info = c_linalg_lapack_dorgqr (_q, tau);
		*q = _q;
	}

	return info;
}

int
c_linalg_QR_unpack (const c_matrix *qr, const c_vector *tau, c_matrix **q, c_matrix **r, bool econ)
{
	int		info;
	info = (econ) ? _c_linalg_QR_unpack_econ (qr, tau, q, r) : _c_linalg_QR_unpack (qr, tau, q, r);
	return info;
}

int
c_linalg_QR_solve (c_matrix *a, c_vector *b)
{
	int			info;
	c_matrix	*x;

	if (c_matrix_is_empty (a)) c_error ("c_linalg_QR_solve", "matrix is empty.");
	if (c_vector_is_empty (b)) c_error ("c_linalg_QR_solve", "vector is empty.");
	if (b->stride != 1) c_error ("c_linalg_QR_solve", "cannot tread vector with stride.");
	if (a->size1 != b->size) c_error ("c_linalg_QR_solve", "vector and matrix size dose not match.");
	if (b->size < a->size2) c_vector_realloc (a->size2, b, b->size);

	x = c_matrix_view_array (b->size, 1, b->size, b->data);
	info = c_linalg_lapack_dgels ('N', a, x);
	c_matrix_free (x);

	if (info == 0 && b->size != a->size2) b->size = a->size2;

	return info;
}

int
c_linalg_lsQ_solve (double rcond, c_matrix *a, c_vector *b, c_vector_int **p, int *rank)
{
	int				info;
	int				_rank;
	c_vector_int	*_p = NULL;
	c_matrix		*x;

	if (c_vector_is_empty (b)) c_error ("c_linalg_lsQ_solve", "vector is empty.");
	if (c_matrix_is_empty (a)) c_error ("c_linalg_lsQ_solve", "matrix is empty.");
	if (b->stride != 1) c_error ("c_linalg_lsQ_solve", "cannot tread vector with stride.");
	if (a->size1 != b->size) c_error ("c_linalg_lsQ_solve", "vector and matrix size dose not match.");
	if (b->size < a->size2) c_vector_realloc (a->size2, b, b->size);

	x = c_matrix_view_array (b->size, 1, b->size, b->data);
	info = c_linalg_lapack_dgelsy (rcond, a, x, &_p, &_rank);
	c_matrix_free (x);

	if (info == 0 && b->size != a->size2) b->size = a->size2;

	if (p) *p = _p;
	else if (!c_vector_int_is_empty (_p)) c_vector_int_free (_p);

	if (rank) *rank = _rank;

	return info;
}

/* solve R *x = Q^T * y */
void
c_linalg_QR_Rsolve (c_matrix *r, c_vector *qty)
{
	char	uplo;
	char	trans;
	char	diag;
	int		n;
	int		lda;
	int		incy;

	if (c_matrix_is_empty (r)) c_error ("c_linalg_QR_Rsolve", "matrix is empty.");
	if (c_vector_is_empty (qty)) c_error ("c_linalg_QR_Rsolve", "vector is empty.");
	if (qty->stride != 1) c_error ("c_linalg_QR_Rsolve", "cannot tread vector with stride.");
	if (qty->size != r->size1) c_error ("c_linalg_QR_Rsolve", "vector and matrix size dose not match.");
	if (r->size2 > qty->size) {
		qty->data = (double *) realloc (qty->data, r->size2 * sizeof (double));
		qty->tsize = r->size2 * qty->stride;
	}

	uplo = 'U';
	trans = 'N';
	diag = 'N';
	n = (int) r->size1;
	lda = (int) r->lda;
	incy = (int) qty->stride;
	dtrsv_ (&uplo, &trans, &diag, &n, r->data, &lda, qty->data, &incy);
	if (qty->size != r->size2) qty->size = r->size2;
	return;
}

/* solve R^T *(Q^T * x) = y */
void
c_linalg_QR_RTsolve (c_matrix *r, c_vector *y)
{
	char	uplo;
	char	trans;
	char	diag;
	int		n;
	int		lda;
	int		incy;

	if (c_matrix_is_empty (r)) c_error ("c_linalg_QR_Rsolve", "matrix is empty.");
	if (c_vector_is_empty (y)) c_error ("c_linalg_QR_Rsolve", "vector is empty.");
	if (y->stride != 1) c_error ("c_linalg_QR_RTsolve", "cannot tread vector with stride.");
	if (y->size != r->size2) c_error ("c_linalg_QR_Rsolve", "vector and matrix size dose not match.");
	if (r->size1 > y->size) {
		y->data = (double *) realloc (y->data, r->size1 * sizeof (double));
		y->tsize = r->size1 * y->stride;
	}

	uplo = 'U';
	trans = 'T';
	diag = 'N';
	n = (int) r->size1;
	lda = (int) r->lda;
	incy = (int) y->stride;
	dtrsv_ (&uplo, &trans, &diag, &n, r->data, &lda, y->data, &incy);
	if (y->size != r->size1) y->size = r->size1;
	return;
}

/*** qr 1-rank update a + u * v' : u and v are vector ***/
void
c_linalg_QR_1up (c_matrix *q, c_matrix *r, const c_vector *u, const c_vector *v)
{
	int		m, n, k, ldq, ldr;
	double	*w;

	if (c_matrix_is_empty (q)) c_error ("c_linalg_QR_1up", "matrix *q is empty.");
	if (c_matrix_is_empty (r)) c_error ("c_linalg_QR_1up", "matrix *r is empty.");
	if (c_vector_is_empty (u)) c_error ("c_linalg_QR_1up", "vector *u is empty.");
	if (c_vector_is_empty (v)) c_error ("c_linalg_QR_1up", "vector *v is empty.");
	if (q->size2 != r->size1) c_error ("c_linalg_QR_1up", "matrix size dose not match.");
	if (u->stride != 1 || v->stride != 1) c_error ("c_linalg_QR_1up", "cannot tread vector with stride.");

	m = q->size1;
	n = r->size2;
	if (c_matrix_is_square (q)) k = q->size2;
	else k = C_MIN (m, n);
	w = (double *) malloc (2 * k * sizeof (double));

	ldq = q->lda;
	ldr = r->lda;

	dqr1up_ (&m, &n, &k, q->data, &ldq, r->data, &ldr, u->data, v->data, w);
	free (w);

	return;
}

/*** insert ***/
void
c_linalg_QR_colinsert (c_matrix *q, c_matrix *r, const size_t index, const c_vector *u)
{
	int			j, m, n, k, ldq, ldr;
	double		*w;

	if (c_matrix_is_empty (q)) c_error ("c_linalg_QR_colinsert", "matrix *q is empty.");
	if (c_matrix_is_empty (r)) c_error ("c_linalg_QR_colinsert", "matrix *r is empty.");
	if (c_vector_is_empty (u)) c_error ("c_linalg_QR_colinsert", "vector *u is empty.");
	if (u->size != q->size1) c_error ("c_linalg_QR_colinsert", "vector and matrix size dose not match..");
	if (q->size2 != r->size1) c_error ("c_linalg_QR_colinsert", "matrix size dose not match..");
	if (index < 0 || r->size2 < index) c_error ("c_linalg_QR_colinsert", "index out of range.");
	if (u->stride != 1) c_error ("c_linalg_QR_colinsert", "cannot tread vector with stride.");

	m = q->size1;
	n = r->size2;
	k = C_MIN (m, n);
	ldq = q->lda;
	if (q->size1 > q->size2) {	// economy mode
		/*=
		 *= | Q11 Q12 |    | Q11 Q12 D0 |
		 *= | Q21 Q22 | -> | Q21 Q22 D0 |
		 *= | Q31 Q32 |    | Q31 Q32 D0 |
		 *=
		 *= | R11 R12 |    | R11 R12 D0 |
		 *= | D0  R22 | -> | D0  R22 D0 |
		 *=                | D0  D0  D0 |
		 */
		c_matrix_add_rowcols (q, 0, 1);
		c_matrix_add_rowcols (r, 1, 1);
	} else c_matrix_add_rowcols (r, 0, 1);

	ldr = r->lda;

	j = index + 1;	// differ from fortran to C

	w = (double *) malloc (k * sizeof (double));
	dqrinc_ (&m, &n, &k, q->data, &ldq, r->data, &ldr, &j, u->data, w);
	free (w);

	return;
}

void
c_linalg_QR_rowinsert (c_matrix *q, c_matrix *r, const size_t index, const c_vector *u)
{
	int			j, m, n, k, ldq, ldr;
	double		*w;

	if (c_matrix_is_empty (q)) c_error ("c_linalg_QR_rowinsert", "matrix *q is empty.");
	if (c_matrix_is_empty (r)) c_error ("c_linalg_QR_rowinsert", "matrix *r is empty.");
	if (c_vector_is_empty (u)) c_error ("c_linalg_QR_rowinsert", "vector *u is empty.");
	if (u->size != r->size2) c_error ("c_linalg_QR_rowinsert", "matrix size dose not match..");
	if (q->size2 != r->size1) c_error ("c_linalg_QR_rowinsert", "matrix size dose not match..");
	if (index < 0 || q->size1 < index) c_error ("c_linalg_QR_rowinsert", "index out of range.");
	if (!c_matrix_is_square (q) && c_matrix_is_square (r))
		c_error ("c_linalg_QR_rowinsert", "rowinsert cannot treat QR economy mode.");
	if (u->stride != 1) c_error ("c_linalg_QR_rowinsert", "cannot tread vector with stride.");

	m = q->size1;
	n = r->size2;

	/*=
	 *= | Q11 Q12 |        | Q11 Q12 D0 |
	 *= | Q21 Q22 |     -> | Q21 Q22 D0 |
	 *=                    | D0  D0  D0 |
	 *=
	 *= | R11 R12 R13 |    | R11 R12 R13 |
	 *= | D0  R22 R23 | -> | D0  R22 R23 |
	 *=                    | D0  D0  D0  |
	 */
	c_matrix_add_rowcols (q, 1, 1);
	c_matrix_add_rowcols (r, 1, 0);

	ldq = q->lda;
	ldr = r->lda;

	j = index + 1;	// differ from fortran to C

	k = C_MIN (m, n);
	w = (double *) malloc (k * sizeof (double));
	dqrinr_ (&m, &n, q->data, &ldq, r->data, &ldr, &j, u->data, w);
	free (w);

	return;
}

/*** delete ***/
void
c_linalg_QR_coldelete (c_matrix *q, c_matrix *r, const size_t index)
{
	int			j, m, n, k, ldq, ldr;
	double		*w;

	if (c_matrix_is_empty (q)) c_error ("c_linalg_QR_coldelete", "matrix *q is empty.");
	if (c_matrix_is_empty (r)) c_error ("c_linalg_QR_coldelete", "matrix *r is empty.");
	if (q->size2 != r->size1) c_error ("c_linalg_QR_coldelete", "matrix size dose not match..");
	if (index < 0 || r->size2 <= index) c_error ("c_linalg_QR_coldelete", "index out of range.");

	m = q->size1;
	n = r->size2;
	k = C_MIN (m, n);
	ldq = q->lda;
	ldr = k;

	j = index + 1;	// differ from fortran to C

	w = (double *) malloc ((k - 1) * sizeof (double));
	dqrdec_ (&m, &n, &k, q->data, &ldq, r->data, &ldr, &j, w);
	free (w);

	if (q->size1 >= q->size2) {
		/*=
		 *= | Q11 Q12 xx |    | Q11 Q12 |
		 *= | Q21 Q22 xx | -> | Q21 Q22 |
		 *= | Q31 Q32 xx |    | Q31 Q32 |
		 *=
		 *= | R11 R12 xx |    | R11 R12 |
		 *= | D0  R22 xx | -> | D0  R22 |
		 *= | D0  D0  xx |
		 */
		c_matrix_remove_rowcols (q, 0, 1);
		c_matrix_remove_rowcols (r, 1, 1);
	} else c_matrix_remove_rowcols (r, 0, 1);

	return;
}

void
c_linalg_QR_rowdelete (c_matrix *q, c_matrix *r, const size_t index)
{
	int			j, m, n, ldq, ldr;
	double		*w;

	if (c_matrix_is_empty (q)) c_error ("c_linalg_QR_rowdelete", "matrix *q is empty.");
	if (c_matrix_is_empty (r)) c_error ("c_linalg_QR_rowdelete", "matrix *r is empty.");
	if (q->size2 != r->size1) c_error ("c_linalg_QR_rowdelete", "matrix size dose not match..");
	if (index < 0 || q->size1 <= index) c_error ("c_linalg_QR_rowdelete", "index out of range.");
	if (!c_matrix_is_square (q) && c_matrix_is_square (r))
		c_error ("c_linalg_QR_rowdelete", "rowdelete cannot treat QR economy mode.");

	m = q->size1;
	n = r->size2;
	ldq = q->lda;
	ldr = r->lda;

	j = index + 1;	// differ from fortran to C

	w = (double *) malloc (2 * m * sizeof (double));
	dqrder_ (&m, &n, q->data, &ldq, r->data, &ldr, &j, w);
	free (w);

	/*=
	 *= | Q11 Q12 xx  |     | Q11 Q12 |
	 *= | Q21 Q22 xx  |  -> | Q21 Q22 |
	 *= | xx  xx  xx  |
	 *=
	 *= | R11 R12 R13 |     | R11 R12 R13 |
	 *= | D0  R22 R23 |  -> | D0  R22 R23 |
	 *= | D0  D0  xx  |
	 *=
	 */
	c_matrix_remove_rowcols (q, 1, 1);
	c_matrix_remove_rowcols (r, 1, 0);

	return;
}
