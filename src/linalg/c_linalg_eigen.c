#include <clinalg.h>

#include "private.h"

/* eigen values of real symmetric matrix */
int
c_linalg_eigen_symm (c_matrix *a, c_vector **e)
{
	int			info;
	c_vector	*_e;
	c_vector	**__e;

	if (c_matrix_is_empty (a)) c_error ("c_linalg_eigen_symm", "matrix a is empty");
	if (!c_matrix_is_square (a)) c_error ("c_linalg_eigen_symm", "matrix a must be square");

	if (e) __e = &_e;
	else __e = NULL;
	info = c_linalg_lapack_dsyev ('N', 'U', a, __e);

	if (e) *e = _e;
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
c_linalg_eigen_general (c_matrix *a, c_vector **er, c_vector **ei, c_matrix **vl, c_matrix **vr)
{
	int			info;
	c_vector	*_er;
	c_vector	*_ei;
	c_matrix	*_vl;
	c_matrix	*_vr;

	c_vector	**__er;
	c_vector	**__ei;
	c_matrix	**__vl;
	c_matrix	**__vr;

	char		jobvl = 'V';
	char		jobvr = 'V';

	if (c_matrix_is_empty (a)) c_error ("c_linalg_eigen_general", "matrix a is empty");

	if (er) __er = &_er;
	else __er = NULL;

	if (ei) __ei = &_ei;
	else __ei = NULL;

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

	info = c_linalg_lapack_dgeev (jobvl, jobvr, a, __er, __ei, __vl, __vr);

	if (er) *er = _er;
	if (ei) *ei = _ei;
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
