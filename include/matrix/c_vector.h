/*
 * c_vector.h
 *
 *  Created on: 2014/04/03
 *      Author: utsugi
 */

#ifndef C_VECTOR_H_
#define C_VECTOR_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#ifndef C_MAX
#define C_MAX(a, b)	((a) >= (b)) ? (a) : (b)
#endif

#ifndef C_MIN
#define C_MIN(a, b)	((a) >= (b)) ? (b) : (a)
#endif

#ifndef INDEX_OF_VECTOR
#define INDEX_OF_VECTOR(v, i) ((i) * (v->stride))
#endif

typedef struct s_c_vector	c_vector;

struct s_c_vector {
	size_t		size;
	size_t		stride;
	size_t		tsize;
	bool		owner;
	double		*data;
};

typedef struct s_c_vector_int	c_vector_int;

struct s_c_vector_int {
	size_t		size;
	size_t		stride;
	size_t		tsize;
	bool		owner;
	int			*data;
};

c_vector		*c_vector_alloc (const size_t size);
void			c_vector_realloc (c_vector *x, const size_t size);
c_vector		*c_vector_view_array (const size_t size, const size_t stride, double *data);
bool			c_vector_is_empty (const c_vector *x);
void			c_vector_free (c_vector *x);
void			c_vector_set (c_vector *x, const int i, const double val);
double			c_vector_get (const c_vector *x, const int i);

void			c_vector_memcpy (c_vector *dest, const c_vector *src);
void			c_vector_set_all (c_vector *x, const double val);
void			c_vector_set_zero (c_vector *x);

c_vector		*c_vector_subvector (const size_t size, const c_vector *x);

void			c_vector_fprintf (FILE *stream, const c_vector *x, const char *format);

c_vector_int	*c_vector_int_alloc (const size_t size);
c_vector_int	*c_vector_int_view_array (const size_t size, const size_t stride, int *data);
bool			c_vector_int_is_empty (const c_vector_int *v);
void			c_vector_int_free (c_vector_int *v);
void			c_vector_int_set (c_vector_int *v, const int i, int val);
int				c_vector_int_get (const c_vector_int *v, const int i);

void			c_vector_int_memcpy (c_vector_int *dest, const c_vector_int *src);

void			c_vector_int_set_zero (c_vector_int *v);
void			c_vector_int_fprintf (FILE *stream, const c_vector_int *v, const char *format);

#ifdef __cplusplus
}
#endif

#endif /* C_VECTOR_H_ */
