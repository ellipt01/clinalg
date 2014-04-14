/*
 * c_linalg_qr.c
 *
 *  Created on: 2014/04/11
 *      Author: utsugi
 */

#include <c_linalg.h>

/* c_linalg_util.c */
extern void	c_error (const char * function_name, const char *error_msg);

/* blas */
extern void	dcopy_ (long *n, double *x, long *incx, double *y, long *incy);

/* lapack */
extern void	dgeqrf_ (long *m, long *n, double *data, long *lda, double *tau, double *work, long *lwork, long *info);
extern void	dgeqp3_ (long *m, long *n, double *a, long *lda, long *jpvt, double *tau, double *work, long *lwork, long *info);
extern void	dorgqr_ (long *m, long *n, long *k, double *data, long *lda, double *tau, double *work, long *lwork, long *info);
extern void	dgels_  (char *trans, long *m, long *n, long *nrhs, double *a_data, long *lda, double *b_data, long *ldb, double *w, long *lwork, long *info);
extern void	dgelsy_ (long *m, long *n, long *nrhs, double *a, long *lda, double *b, long *ldb, long *jpvt, double *rcond, long *rank, double *work, long *lwork, long *info);
extern void	dtrsv_ (char *uplo, char *trans, char *diag, long *n, double *r, long *lda, double *y, long *incy);

/* qrupdate*/
extern void	dqr1up_ (int *m, int *n, int *k, double *Q, int *ldq, double *R, int *ldr, double *u, double *v, double *w);
extern void	dqrinc_ (int *m, int *n, int *k, double *Q, int *ldq, double *R, int *ldr, int *j, double *x, double *w);
extern void	dqrinr_ (int *m, int *n, double *Q, int *ldq, double *R, int *ldr, int *j, double *x, double *w);
extern void	dqrdec_ (int *m, int *n, int *k, double *Q, int *ldq, double *R, int *ldr, int *j, double *w);
extern void	dqrder_ (int *m, int *n, double *Q, int *ldq, double *R, int *ldr, int *j, double *w);

int
c_linalg_lapack_dgeqrf (c_matrix *a, c_vector **tau)
{
	long   	info;
	long		m;
	long		n;
	long		lda;

	size_t		min_mn;
	size_t		ltau;

	double		wkopt;
	long		lwork;
	double		*work;

	c_vector		*_tau;

	if (c_matrix_is_empty (a)) c_error ("c_linalg_lapack_dgeqrf", "matrix is empty.");

	m = (long) a->size1;
	n = (long) a->size2;
	lda = (long) a->lda;

	min_mn = (size_t) C_MIN (a->size1, a->size2);
	ltau = (size_t) C_MAX (min_mn, 1);

	_tau = c_vector_alloc (ltau);

 	lwork = -1;
 	dgeqrf_ (&m, &n, a->data, &lda, _tau->data, &wkopt, &lwork, &info);

 	lwork = (long) wkopt;
	if ((int) info != 0 || lwork <= 0) c_error ("c_linalg_lapack_dgeqrf", "failed to query size of workspace.");
	if ((work = (double *) malloc (lwork * sizeof (double))) == NULL)
		c_error ("c_linalg_lapack_dgeqrf", "cannot allocate memory for workspace.");
	dgeqrf_ (&m, &n, a->data, &lda, _tau->data, work, &lwork, &info);
	free (work);

	if (tau) *tau = _tau;
	else c_vector_free (_tau);

	return (int) info;
}

int
c_linalg_lapack_dgeqp3 (c_matrix *a, c_vector **tau, long **p)
{
	long   	info;
	long		m;
	long		n;
	long		lda;

	size_t		min_mn;
	size_t		ltau;

	double		wkopt;
	long		lwork;
	double		*work;

	c_vector	*_tau;
	long		*_p;

	if (c_matrix_is_empty (a)) c_error ("c_linalg_lapack_dgeqp3", "matrix is empty.");

	m = (long) a->size1;
	n = (long) a->size2;
	lda = (long) a->lda;

	min_mn = (size_t) C_MIN (a->size1, a->size2);
	ltau = (size_t) C_MAX (min_mn, 1);

	_tau = c_vector_alloc (ltau);
	_p = (long *) malloc (min_mn * sizeof (long));

	lwork = -1;
	dgeqp3_ (&m, &n, a->data, &lda, _p, _tau->data, &wkopt, &lwork, &info);

	lwork = (int) wkopt;
	if ((int) info != 0 || lwork <= 0) c_error ("c_linalg_lapack_dgeqp3", "failed to query size of workspace.");
	if ((work = (double *) malloc (lwork * sizeof (double))) == NULL)
		c_error ("c_linalg_lapack_dgeqp3", "cannot allocate memory for workspace.");
	dgeqp3_ (&m, &n, a->data, &lda, _p, _tau->data, work, &lwork, &info);
	free (work);

	if (p) *p = _p;
	else free (_p);

	if (tau) *tau = _tau;
	else c_vector_free (_tau);

	return (int) info;
}

int
c_linalg_lapack_dorgqr (c_matrix *qr, c_vector *tau)
{
	long   	info;
 	long		m;
 	long		n;
 	long		min_mn;
 	long		k;
	long		lda;
	long		lwork;
	double		wkopt;
	double		*work;

	if (c_matrix_is_empty (qr)) c_error ("c_linalg_lapack_dorgqr", "matrix is empty.");
	if (c_vector_is_empty (tau)) c_error ("c_linalg_lapack_dorgqr", "vector *tau is empty.");

	m = (long) qr->size1;
	n = (long) qr->size2;
	min_mn = (long) C_MIN (m, n);
	k = (long) C_MAX (min_mn, 1);
	lda = (long) qr->lda;

 	if (tau->size != k) c_error ("c_linalg_lapack_dorgqr", "tau->size must be equal to MAX (MIN (m, n), 1).");

 	lwork = -1;
	dorgqr_ (&m, &min_mn, &k, qr->data, &lda, tau->data, &wkopt, &lwork, &info);
	lwork = (long) wkopt;
	if ((int) info != 0 || lwork <= 0) c_error ("c_linalg_lapack_dorgqr", "failed to query size of workspace.");
	if ((work = (double *) malloc (lwork * sizeof (double))) == NULL)
		c_error ("c_linalg_lapack_dorgqr", "cannot allocate memory for workspace.");
	dorgqr_ (&m, &min_mn, &k, qr->data, &lda, tau->data, work, &lwork, &info);
	free (work);

	if (qr->size1 < qr->size2) qr->size2 = qr->size1;

	return (int) info;
}

int
c_linalg_lapack_dgels (char trans, c_matrix *qr, c_matrix *b)
{
	long		info;
	long		m;
	long		n;
	long		nrhs;
	long		lda;
	long		ldb;
	double		wkopt;
	double		*work;
	long		lwork;

	if (c_matrix_is_empty(qr)) c_error ("c_linalg_lapack_dgels", "input matrix *qr is empty.");
	if (c_matrix_is_empty(b)) c_error ("c_linalg_lapack_dgels", "input matrix *b is empty.");

	m = (long) qr->size1;
	n = (long) qr->size2;
	nrhs = (long) b->size2;
	lda = (long) qr->lda;
	ldb = (long) C_MAX (1, C_MAX (m, n));

	lwork = -1;
	dgels_ (&trans, &m, &n, &nrhs, qr->data, &lda, b->data, &ldb, &wkopt, &lwork, &info);
	lwork = (long) wkopt;
	if ((int) info != 0 || lwork <= 0) c_error ("c_linalg_lapack_dgels", "failed to query size of workspace.");
	if ((work = (double *) malloc (lwork * sizeof (double))) == NULL)
		c_error ("c_linalg_lapack_dgels", "cannot allocate memory for workspace.");
	dgels_ (&trans, &m, &n, &nrhs, qr->data, &lda, b->data, &ldb, work, &lwork, &info);
	free (work);

	return (int) info;
}

int
c_linalg_lapack_dgelsy (double rcond, c_matrix *qr, c_matrix *b, long **p, int *rank)
{
	long		info;
	long		m;
	long		n;
	long		nrhs;
	long		lda;
	long		ldb;
	long		lrank;
	double		wkopt;
	double		*work;
	long		lwork;
	long		*_p;

	if (c_matrix_is_empty (qr)) c_error ("c_linalg_lapack_dgelsy", "input matrix *qr is empty");
	if (c_matrix_is_empty (b)) c_error ("c_linalg_lapack_dgelsy", "input matrix *b is empty");
	if (qr->size1 < qr->size2 && b->size1 < qr->size2) {
		double	*data = (double *) realloc (b->data, qr->size2 * b->size2);
		if (!data) c_error ("c_linalg_lapack_dgelsy", "cannot reallocate matrix.");
		b->data = data;
	}

	m = (long) qr->size1;
	n = (long) qr->size2;
	nrhs = (long) b->size2;
	lda = (long) qr->lda;
	ldb = (long) C_MAX (1, C_MAX (m, n));
	_p = (long *) malloc (qr->size2 * sizeof (long));

	lwork = -1;
	dgelsy_ (&m, &n, &nrhs, qr->data, &lda, b->data, &ldb, _p, &rcond, &lrank, &wkopt, &lwork, &info);

	lwork = (long) wkopt;
	if ((int) info != 0 || lwork <= 0) c_error ("c_linalg_lapack_dgelsy", "failed to query workspace");
	if ((work = (double *) malloc (lwork * sizeof (double))) == NULL)
		c_error ("c_linalg_lapack_dgelsy", "cannot allocate memory work");
	dgelsy_ (&m, &n, &nrhs, qr->data, &lda, b->data, &ldb, _p, &rcond, &lrank, work, &lwork, &info);
	free (work);

	if (p) *p = _p;
	else free (_p);

	if (rank) *rank = (int) lrank;

	return (int) info;
}

int
c_linalg_QR_decomp (c_matrix *a, long **p, c_vector **tau)
{
	int			info;
 	long		*_p;
	c_vector	*_tau;

	if (c_matrix_is_empty (a)) c_error ("c_linalg_QR_decomp", "matrix is empty.");

	if (p == NULL) info = c_linalg_lapack_dgeqrf (a, &_tau);
	else info = c_linalg_lapack_dgeqp3 (a, &_tau, &_p);

	if (p) *p = _p;
	if (tau) *tau = _tau;
	else c_vector_free (_tau);

	return info;
}

int
c_linalg_QR_unpack (c_matrix *qr, c_vector *tau)
{
	int		info;

	if (c_matrix_is_empty (qr)) c_error ("c_linalg_QR_unpack", "matrix is empty.");

	info = c_linalg_lapack_dorgqr (qr, tau);

	return info;
}

int
c_linalg_QR_solve (c_matrix *a, c_vector *b)
{
	int			info;
	c_matrix	*x;

	if (c_matrix_is_empty (a)) c_error ("c_linalg_QR_solve", "matrix is empty.");
	if (c_vector_is_empty (b)) c_error ("c_linalg_QR_solve", "vector is empty.");
	if (a->size1 != b->size) c_error ("c_linalg_QR_solve", "vector and matrix size dose not match.");
	if (a->size1 < a->size2 && b->size < a->size2) {
		size_t	size = b->size;
		c_vector_realloc (b, a->size2);
		b->size = size;
	}

	x = c_matrix_view_array (b->size, 1, b->size, b->data);
	info = c_linalg_lapack_dgels ('N', a, x);
	c_matrix_free (x);

	if (info == 0) b->size = a->size2;

	return (int) info;
}

int
c_linalg_lsQ_solve (double rcond, c_matrix *a, c_vector *b, long **p, int *rank)
{
	int			info;
	int			_rank;
	long		*_p;
	c_matrix	*x;

	if (c_vector_is_empty (b)) c_error ("c_linalg_lsQ_solve", "vector is empty.");
	if (c_matrix_is_empty (a)) c_error ("c_linalg_lsQ_solve", "matrix is empty.");
	if (a->size1 != b->size) c_error ("c_linalg_lsQ_solve", "vector and matrix size dose not match.");
	if (a->size1 < a->size2 && b->size < a->size2) {
		size_t	size = b->size;
		c_vector_realloc (b, a->size2);
		b->size = size;
	}

	x = c_matrix_view_array (b->size, 1, b->size, b->data);
	info = c_linalg_lapack_dgelsy (rcond, a, x, &_p, &_rank);
	c_matrix_free (x);

	if (info == 0) b->size = a->size2;

	if (p) *p = _p;
	if (rank) *rank = _rank;

	return info;
}

void
c_linalg_QR_Rsolve (c_matrix *r, c_vector *qty)
{
	char	uplo;
	char	trans;
	char	diag;
	long	n;
	long	lda;
	long	incy;

	if (c_matrix_is_empty (r)) c_error ("c_linalg_QR_Rsolve", "matrix is empty.");
	if (c_vector_is_empty (qty)) c_error ("c_linalg_QR_Rsolve", "vector is empty.");
	if (!c_matrix_is_square (r)) c_error ("c_linalg_QR_Rsolve", "matrix must be square.");
	if (qty->size != r->size1) c_error ("c_linalg_QR_Rsolve", "vector and matrix size dose not match.");

	uplo = 'U';
	trans = 'N';
	diag = 'N';
	n = (long) r->size1;
	lda = (long) r->lda;
	incy = (long) qty->stride;
	dtrsv_ (&uplo, &trans, &diag, &n, r->data, &lda, qty->data, &incy);
	if (qty->size != r->size2) qty->size = r->size2;
	return;
}

/*** qr 1-rank update a + u * v' : u and v are vector ***/
void
c_linalg_QR_1up (c_matrix *q, c_matrix *r, c_vector *u, c_vector *v)
{
	int		m, n, k, ldq, ldr;
	double	*w;

	if (c_matrix_is_empty (q)) c_error ("c_linalg_QR_1up", "matrix *q is empty.");
	if (c_matrix_is_empty (r)) c_error ("c_linalg_QR_1up", "matrix *r is empty.");
	if (c_vector_is_empty (u)) c_error ("c_linalg_QR_1up", "vector *u is empty.");
	if (c_vector_is_empty (v)) c_error ("c_linalg_QR_1up", "vector *v is empty.");
	if (q->size2 != r->size1) c_error ("c_linalg_QR_1up", "matrix size dose not match.");

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
c_linalg_QR_colinsert (c_matrix *q, c_matrix *r, const size_t index, c_vector *u)
{
	int			j, m, n, k, ldq, ldr;
	double		*w;

	if (c_matrix_is_empty (q)) c_error ("c_linalg_QR_colinsert", "matrix *q is empty.");
	if (c_matrix_is_empty (r)) c_error ("c_linalg_QR_colinsert", "matrix *r is empty.");
	if (c_vector_is_empty (u)) c_error ("c_linalg_QR_colinsert", "vector *u is empty.");
	if (q->size2 != r->size1) c_error ("c_linalg_QR_colinsert", "matrix size dose not match..");
	if (index < 0 || r->size2 < index) c_error ("c_linalg_QR_colinsert", "index out of range.");

	m = q->size1;
	n = r->size2;
	k = q->size2;
	ldq = q->lda;
	if (!c_matrix_is_square (q)) {
	/*
		| Q11 Q12 |    | Q11 Q12 D0 |
		| Q21 Q22 | -> | Q21 Q22 D0 |
		| Q31 Q32 |    | Q31 Q32 D0 |

		| R11 R12 |    | R11 R12 D0 |
		| D0  R22 | -> | D0  R22 D0 |
                      | D0  D0  D0 |
	*/
		c_matrix_add_col (q);
		c_matrix_add_row (r);
	}
	c_matrix_add_col (r);
	ldr = r->lda;

	j = index + 1;	// differ from fortran to C

	w = (double *) malloc (k * sizeof (double));
	dqrinc_ (&m, &n, &k, q->data, &ldq, r->data, &ldr, &j, u->data, w);
	free (w);

	return;
}

void
c_linalg_QR_rowinsert (c_matrix *q, c_matrix *r, const size_t index, c_vector *u)
{
	int			j, m, n, k, ldq, ldr;
	double		*w;

	if (c_matrix_is_empty (q)) c_error ("c_linalg_QR_rowinsert", "matrix *q is empty.");
	if (c_matrix_is_empty (r)) c_error ("c_linalg_QR_rowinsert", "matrix *r is empty.");
	if (c_vector_is_empty (u)) c_error ("c_linalg_QR_rowinsert", "vector *u is empty.");
	if (q->size2 != r->size1) c_error ("c_linalg_QR_rowinsert", "matrix size dose not match..");
	if (index < 0 || r->size1 < index) c_error ("c_linalg_QR_rowinsert", "index out of range.");
	if (!c_matrix_is_square (q) && c_matrix_is_square (r))
		c_error ("c_linalg_QR_rowinsert", "rowinsert cannot treat QR economy mode.");

	m = q->size1;
	n = r->size2;
	if (c_matrix_is_square (q)) {
	/*
		| Q11 Q12 |        | Q11 Q12 D0 |
		| Q21 Q22 |     -> | Q21 Q22 D0 |
		                   | D0  D0  D0 |

		| R11 R12 R13 |    | R11 R12 R13 |
		| D0  R22 R23 | -> | D0  R22 R23 |
		                   | D0  D0  D0  |
	*/
		c_matrix_add_col (q);
		c_matrix_add_row (r);
	}
	c_matrix_add_row (q);
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
	k = q->size2;
	ldq = q->lda;
	ldr = r->lda;

	j = index + 1;	// differ from fortran to C

	w = (double *) malloc ((k - 1) * sizeof (double));
	dqrdec_ (&m, &n, &k, q->data, &ldq, r->data, &ldr, &j, w);
	free (w);

	if (!c_matrix_is_square (q)) {
	/*
		| Q11 Q12 xx |    | Q11 Q12 |
		| Q21 Q22 xx | -> | Q21 Q22 |
		| xx  xx  xx |

		| R11 R12 xx |    | R11 R12 |
		| D0  R22 xx | -> | D0  R22 |
       | D0  D0  xx |
	*/
		c_matrix_remove_col (q);
		c_matrix_remove_row (r);
	}
	c_matrix_remove_col (r);

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
	if (index < 0 || r->size1 <= index) c_error ("c_linalg_QR_rowdelete", "index out of range.");
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

	if (c_matrix_is_square (q)) {
	/*
		| Q11 Q12 xx |     | Q11 Q12 |
		| Q21 Q22 xx |  -> | Q21 Q22 |
		| xx  xx  xx |

		| R11 R12 R13 |    | R11 R12 R13 |
		| D0  R22 R23 | -> | D0  R22 R23 |
		| xx  xx  xx  |
	*/
		c_matrix_remove_col (q);
		c_matrix_remove_row (q);
	}
	c_matrix_remove_row (r);

	return;
}
