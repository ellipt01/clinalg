/*
 * c_matrix_double.h
 *
 *  Created on: 2014/04/03
 *      Author: utsugi
 */

#ifndef C_MATRIX_DOUBLE_H_
#define C_MATRIX_DOUBLE_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <c_vector_double.h>

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

void			c_matrix_add_col (c_matrix *a);
void			c_matrix_add_row (c_matrix *a);
void			c_matrix_remove_col (c_matrix *a);
void			c_matrix_remove_row (c_matrix *a);

void			c_matrix_memcpy (c_matrix *dest, const c_matrix *src);
void			c_matrix_set_zero (c_matrix *a);
void			c_matrix_fprintf (FILE *stream, const c_matrix *a, const char *format);

#ifdef __cplusplus
}
#endif

#endif /* C_MATRIX_DOUBLE_H_ */
