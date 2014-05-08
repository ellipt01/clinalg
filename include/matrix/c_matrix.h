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

#ifndef INDEX_OF_MATRIX
#define INDEX_OF_MATRIX(a, i, j) ((i) + (j) * (a->lda))
#endif

#ifndef POINTER_OF_MATRIX
#define POINTER_OF_MATRIX(a, i, j) ((a->data) + (i) + (j) * (a->lda))
#endif

typedef struct s_c_matrix	c_matrix;

struct s_c_matrix {
	size_t		size1;
	size_t		size2;
	size_t		lda;
	size_t		tsize;
	bool		owner;
	double		*data;
};

c_matrix		*c_matrix_alloc (const size_t size1, const size_t size2);
c_matrix		*c_matrix_view_array (const size_t size1, const size_t size2, const size_t lda, double *data);
bool			c_matrix_is_empty (const c_matrix *a);
bool			c_matrix_is_square (const c_matrix *a);
void			c_matrix_free (c_matrix *a);
void			c_matrix_set (c_matrix *a, const int i, const int j, double val);
double			c_matrix_get (const c_matrix *a, const int i, const int j);

void			c_matrix_get_row (c_vector *v, const c_matrix *a, const size_t index);
void			c_matrix_get_col (c_vector *v, const c_matrix *a, const size_t index);
void			c_matrix_set_row (c_matrix *a, const size_t index, const c_vector *v);
void			c_matrix_set_col (c_matrix *a, const size_t index, const c_vector *v);

c_vector		*c_matrix_column (c_matrix *a, int index);
c_vector		*c_matrix_row (c_matrix *a, int index);

void			c_matrix_memcpy (c_matrix *dest, const c_matrix *src);
void			c_matrix_mncopy (c_matrix *dest, const size_t m0, const size_t n0, const size_t m, const size_t n, const c_matrix *src);

void			c_matrix_set_all (c_matrix *a, const double val);
void			c_matrix_set_zero (c_matrix *a);

void			c_matrix_set_diagonal (const c_vector *d, c_matrix *a);
c_vector		*c_matrix_get_diagonal (const c_matrix *a);
c_vector		*c_matrix_diagonal_view_array (const c_matrix *a);

c_matrix		*c_matrix_submatrix (const size_t m0, const size_t n0, const size_t m, const size_t n, const c_matrix *a);

void			c_matrix_fprintf (FILE *stream, const c_matrix *a, const char *format);
void			c_matrix_fprintf2 (FILE *stream, const c_matrix *a, const char *format);
c_matrix		*c_matrix_fread (FILE *stream, const size_t size1, const size_t size2);

#ifdef __cplusplus
}
#endif

#endif /* C_MATRIX_H_ */
