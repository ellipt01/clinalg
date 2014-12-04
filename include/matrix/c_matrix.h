/*
 * c_matrix.h
 *
 *  Created on: 2014/04/03
 *      Author: utsugi
 */

#ifndef C_MATRIX_H_
#define C_MATRIX_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <c_vector.h>

typedef struct s_c_matrix	c_matrix;

struct s_c_matrix {
	int		size1;
	int		size2;
	int		lda;
	int		tsize;
	bool	owner;
	double	*data;
};

c_matrix		*c_matrix_alloc (const int size1, const int size2);
void			c_matrix_realloc (const int tsize, c_matrix *x, const int size1, const int size2);
c_matrix		*c_matrix_view_array (const int size1, const int size2, const int lda, double *data);
bool			c_matrix_is_empty (const c_matrix *a);
bool			c_matrix_is_square (const c_matrix *a);
void			c_matrix_free (c_matrix *a);
void			c_matrix_set (c_matrix *a, const int i, const int j, const double val);
double			c_matrix_get (const c_matrix *a, const int i, const int j);

void			c_matrix_get_row (c_vector *v, const c_matrix *a, const int index);
void			c_matrix_get_col (c_vector *v, const c_matrix *a, const int index);
void			c_matrix_set_row (c_matrix *a, const int index, const c_vector *v);
void			c_matrix_set_col (c_matrix *a, const int index, const c_vector *v);

c_vector		*c_matrix_column (c_matrix *a, const int index);
c_vector		*c_matrix_row (c_matrix *a, const int index);

void			c_matrix_memcpy (c_matrix *dest, const c_matrix *src);
void			c_matrix_mncopy (c_matrix *dest, const int m0, const int n0, const int m, const int n, const c_matrix *src);
void			c_matrix_upper_triangular_memcpy (c_matrix *tr, const c_matrix *a);
void			c_matrix_lower_triangular_memcpy (c_matrix *tr, const c_matrix *a);

void			c_matrix_set_all (c_matrix *a, const double val);
void			c_matrix_set_zero (c_matrix *a);

void			c_matrix_set_diagonal (const c_vector *d, c_matrix *a);
c_vector		*c_matrix_get_diagonal (const c_matrix *a);
c_vector		*c_matrix_diagonal_view_array (const c_matrix *a);

c_matrix		*c_matrix_submatrix (const int m0, const int n0, const int m, const int n, const c_matrix *a);

void			c_matrix_fprintf (FILE *stream, const c_matrix *a, const char *format);
void			c_matrix_fprintf2 (FILE *stream, const c_matrix *a, const char *format);
c_matrix		*c_matrix_fread (FILE *stream, const int size1, const int size2);

#ifdef __cplusplus
}
#endif

#endif /* C_MATRIX_H_ */
