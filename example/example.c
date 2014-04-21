/*
 ============================================================================
 Name        : exampleProgram.c
 Author      : utsugi
 Version     :
 Copyright   : GPL
 Description : Uses shared library to print greeting
               To run the resulting executable the LD_LIBRARY_PATH must be
               set to ${project_loc}/liblarsen/.libs
               Alternatively, libtool creates a wrapper shell script in the
               build directory of this program which can be used to run it.
               Here the script will be called exampleProgram.
 ============================================================================
 */

#include <clinalg.h>

int
main (void)
{
	c_matrix	*a = c_matrix_alloc (10, 10);

	c_matrix_free (a);
	return 0;
}
