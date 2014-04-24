/*
 * c_linalg_lu.c
 *
 *  Created on: 2014/04/03
 *      Author: utsugi
 */

#include <clinalg.h>

/* c_linalg_util.c */
extern void	c_error (const char * function_name, const char *error_msg);

/* lapack */
#ifndef HAVE_LAPACK_H
extern void	dgetrf_ (int *m, int *n, double *data, int *lda, int *ipiv, int *info);
extern void	dgetrs_ (char *trans, int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info);
extern void	dgetri_ (int *n, double *data, int *lda, int *ipiv, double *work, int *lwork, int *info);
#endif
/* qrupdate */
extern void	dlup1up_ (int *m, int *n, double *L, int *ldl, double *R, int *ldr, int *p, double *u, double *v, double *w);

int
c_linalg_lapack_dgetrf (c_matrix *a, c_vector_int **p)
{
	int				info;
 	int				m;
 	int				n;
 	int				lda;
 	size_t			min_mn;
 	c_vector_int	*_p;

 	if (c_matrix_is_empty (a)) c_error ("c_linalg_lapack_dgetrf", "matrix is empty.");

 	m = (int) a->size1;
 	n = (int) a->size2;
 	lda  = (int) a->lda;
 	min_mn = (size_t) C_MIN (a->size1, a->size2);
 	_p = c_vector_int_alloc (min_mn);
	dgetrf_ (&m, &n, a->data, &lda, _p->data, &info);

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

int
c_linalg_LU_decomp (c_matrix *a, c_vector_int **p)
{
	int				info;
	c_vector_int	*_p = NULL;

	if (c_matrix_is_empty (a)) c_error ("c_linalg_LU_decomp", "matrix is empty.");

	info = c_linalg_lapack_dgetrf (a, &_p);

	if (p) *p = _p;
	else if (!c_vector_int_is_empty (_p)) c_vector_int_free (_p);

	return info;
}

void
c_linalg_LU_unpack (const c_matrix *lu, c_matrix **l, c_matrix **u)
{
	int		min_mn;

	if (c_matrix_is_empty (lu)) c_error ("c_linalg_LU_unpack", "matrix is empty.");

	min_mn = (int) C_MIN (lu->size1, lu->size2);

	if (l) {
		int			i;
		c_matrix	*_l = c_matrix_alloc (lu->size1, min_mn);
		c_matrix_set_zero (_l);
		c_matrix_lower_triangular_memcpy (_l, lu);
		for (i = 0; i < min_mn; i++) c_matrix_set (_l, i, i, 1.);
		*l = _l;
	}

	if (u) {
		c_matrix	*_u = c_matrix_alloc (min_mn, lu->size2);
		c_matrix_set_zero (_u);
		c_matrix_upper_triangular_memcpy (_u, lu);
		*u = _u;
	}
	return;
}

int
c_linalg_LU_solve (c_matrix *lu, c_vector *b, c_vector_int *p)
{
	int  		info;
	c_matrix	*c;

	if (c_matrix_is_empty (lu)) c_error ("c_linalg_LU_solve", "matrix is empty.");
	if (c_vector_is_empty (b)) c_error ("c_linalg_LU_solve", "vector is empty.");
	if (!p) c_error ("c_linalg_LU_solve", "permutation is empty.");
	if (b->size != lu->size1) c_error ("c_linalg_LU_solve", "matrix and vector size dose not match.");
	if (!c_matrix_is_square (lu)) c_error ("c_linalg_LU_solve", "matrix must be square.");

	c = c_matrix_view_array (b->size, 1, b->size, b->data);
	info = c_linalg_lapack_dgetrs ('N', lu, c, p);
	c_matrix_free (c);

	return info;
}

int
c_linalg_LU_invert (c_matrix *lu, c_vector_int *p)
{
	int		info;

	if (c_matrix_is_empty (lu)) c_error ("c_linalg_LU_invert", "matrix is empty.");
	if (!c_matrix_is_square (lu)) c_error ("c_linalg_LU_invert", "matrix must be square.");

	info = c_linalg_lapack_dgetri (lu, p);

	return info;
}

void
c_linalg_LU_1up (c_matrix *l, c_matrix *u, c_vector_int *p, c_vector *s, c_vector *t)
{
	int			m;
	int			n;
	int			ldl;
	int			ldu;
	double		*w;

	if (c_matrix_is_empty (l)) c_error ("c_linalg_LU_1up", "matrix is empty.");
	if (c_matrix_is_empty (u)) c_error ("c_linalg_LU_1up", "matrix is empty.");
	if (c_vector_is_empty (s)) c_error ("c_linalg_LU_1up", "vector *s is empty.");
	if (c_vector_is_empty (t)) c_error ("c_linalg_LU_1up", "vector *t is empty.");
	if (c_vector_int_is_empty (p)) c_error ("c_linalg_LU_1up", "permulation is empty.");
	if (s->size != l->size1) c_error ("c_linalg_LU_1up", "vector and matrix size dose not match.");
	if (t->size != u->size2) c_error ("c_linalg_LU_1up", "vector and matrix size dose not match.");

	m = (int) l->size1;
	n = (int) u->size2;
	ldl = (int) l->lda;
	ldu = (int) u->lda;
	w = (double *) malloc (l->size1 * sizeof (double));
	dlup1up_ (&m, &n, l->data, &ldl, u->data, &ldu, p->data, s->data, t->data, w);
	free (w);

	return;
}
