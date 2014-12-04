/*
 * c_linalg_utils.c
 *
 *  Created on: 2014/04/10
 *      Author: utsugi
 */

#include <clinalg.h>

const int		izero = 0;
const int		ione  = 1;
const double	dzero = 0.;
const double	done  = 1.;
const double	dmone = -1.;

void
c_error (const char * function_name, const char *error_msg)
{
	fprintf (stderr, "ERROR: %s: %s\n", function_name, error_msg);
	exit (1);
}
