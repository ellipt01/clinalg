#include <clinalg.h>

#include "private.h"

/* eigen values of real symmetric matrix */
int
c_linalg_eigen_symm (c_matrix *a, c_vector **w)
{
	int			info;
	c_vector	*_w;
	c_vector	**__w;

	if (c_matrix_is_empty (a)) c_error ("c_linalg_eigen_symm", "matrix a is empty");
	if (!c_matrix_is_square (a)) c_error ("c_linalg_eigen_symm", "matrix a must be square");

	if (w) __w = &_w;
	else __w = NULL;
	info = c_linalg_lapack_dsyev ('N', 'U', a, __w);

	if (w) *w = _w;

	return info;
}

/* eigen values of real symmetric matrix */
int
c_linalg_eigen_symm_d (c_matrix *a, c_vector **e)
{
	int			info;
	c_vector	*_e;
	c_vector	**__e;

	if (c_matrix_is_empty (a)) c_error ("c_linalg_eigen_symm_d", "matrix a is empty");
	if (!c_matrix_is_square (a)) c_error ("c_linalg_eigen_symm_d", "matrix a must be square");

	if (e) __e = &_e;
	else __e = NULL;
	info = c_linalg_lapack_dsyevd ('N', 'U', a, __e);

	if (e) *e = _e;
	return info;
}

/* eigen values and eigen vectors of real general matrix */
int
c_linalg_eigen_general (c_matrix *a, c_vector **wr, c_vector **wi, c_matrix **vl, c_matrix **vr)
{
	int			info;
	c_vector	*_wr;
	c_vector	*_wi;
	c_matrix	*_vl;
	c_matrix	*_vr;

	c_vector	**__wr;
	c_vector	**__wi;
	c_matrix	**__vl;
	c_matrix	**__vr;

	char		jobvl = 'V';
	char		jobvr = 'V';

	if (c_matrix_is_empty (a)) c_error ("c_linalg_eigen_general", "matrix a is empty");

	if (wr) __wr = &_wr;
	else __wr = NULL;

	if (wi) __wi = &_wi;
	else __wi = NULL;

	if (vl) __vl = &_vl;
	else {
		jobvl = 'N';
		__vl = NULL;
	}

	if (vr) __vr = &_vr;
	else {
		jobvr = 'N';
		__vr = NULL;
	}

	info = c_linalg_lapack_dgeev (jobvl, jobvr, a, __wr, __wi, __vl, __vr);

	if (wr) *wr = _wr;
	if (wi) *wi = _wi;
	if (vl) *vl = _vl;
	if (vr) *vr = _vr;
	return info;
}

/* schur decomposition */
static bool
_select (const double *ar, const double *ai)
{
	if (*ai == 0.) return true;
	return false;
}

int
c_linalg_schur_decomp (c_matrix *a, c_vector **wr, c_vector **wi, c_matrix **s)
{
	int			info;
	c_vector	*_wr;
	c_vector	*_wi;
	c_matrix	*_s;

	if (c_matrix_is_empty (a)) c_error ("c_linalg_schur_decomp", "matrix a is empty");
	if (!c_matrix_is_square (a)) c_error ("c_linalg_schur_decomp", "matrix a must be square");

	info = c_linalg_lapack_dgees ('V', 'S', a, &_wr, &_wi, &_s, &_select);
	if (wr) *wr = _wr;
	else c_vector_free (_wr);

	if (wi) *wi = _wi;
	else c_vector_free (_wi);

	if (s) *s = _s;
	else c_matrix_free (_s);

	return info;
}
