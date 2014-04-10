/*
 * c_linalg_utils.c
 *
 *  Created on: 2014/04/10
 *      Author: utsugi
 */

#include <c_linalg.h>

void
c_error (const char * function_name, const char *error_msg)
{
	fprintf (stderr, "ERROR: %s: %s\n", function_name, error_msg);
	exit (1);
}

