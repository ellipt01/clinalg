/*
 * c_vector_int.h
 *
 *  Created on: 2014/04/03
 *      Author: utsugi
 */

#ifndef C_VECTOR_INT_H_
#define C_VECTOR_INT_H_

#ifdef __cplusplus
extern "C" {
#endif

typedef struct s_c_vector_int	c_vector_int;

struct s_c_vector_int {
	size_t		size;
	size_t		stride;
	size_t		tsize;
	bool		owner;
	int			*data;
};

c_vector_int		*c_vector_int_alloc (const size_t size);
c_vector_int		*c_vector_int_view_array (const size_t size, const size_t stride, int *data);
bool				c_vector_int_is_empty (const c_vector_int *v);
void				c_vector_int_free (c_vector_int *v);
void				c_vector_int_set (c_vector_int *v, const int i, int val);
int					c_vector_int_get (const c_vector_int *v, const int i);
void				c_vector_int_memcpy (c_vector_int *dest, const c_vector_int *src);
void				c_vector_int_set_zero (c_vector_int *v);
void				c_vector_int_fprintf (FILE *stream, const c_vector_int *v, const char *format);

#ifdef __cplusplus
}
#endif

#endif /* C_VECTOR_INT_H_ */
