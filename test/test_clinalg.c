/*
 * test.c
 *
 *  Created on: 2014/04/14
 *      Author: utsugi
 */

#include <c_linalg.h>
#include "test_clinalg.h"

const int	n_test_cholesky_func = 4;

const char	*test_cholesky_func_name[] = {
	"test_cholesky_decomp",
	"test_cholesky_svx   ",
	"test_cholesky_insert",
	"test_cholesky_delete"
};

bool	(*test_cholesky_func_ptr[]) (void) = {
	test_cholesky_decomp,
	test_cholesky_svx,
	test_cholesky_insert,
	test_cholesky_delete
};

const int	n_test_QR_func = 5;

const char	*test_QR_func_name[] = {
	"test_QR_decomp",
	"test_QR_solve ",
	"test_lsQ_solve ",
	"test_QR_Rsolve",
	"test_QR_delete"
};

bool	(*test_QR_func_ptr[]) (void) = {
	test_QR_decomp,
	test_QR_solve,
	test_lsQ_solve,
	test_QR_Rsolve,
	test_cholesky_delete
};

int
main (void)
{
	int		i;
	bool	success = true;

	fprintf (stderr, "*** test_cholesky ***\n");
	for (i = 0; i < n_test_cholesky_func; i++) {
		bool	status = test_cholesky_func_ptr[i] ();
		fprintf (stderr, "\t %s ...", test_cholesky_func_name[i]);
		fprintf (stderr, status ? "SUCCESS\n" : "FAILED\n");
		if (success && !status) success = false;
	}

	fprintf (stderr, "*** test_QR ***\n");
	for (i = 0; i < n_test_QR_func; i++) {
		bool	status = test_QR_func_ptr[i] ();
		fprintf (stderr, "\t %s ...", test_QR_func_name[i]);
		fprintf (stderr, status ? "SUCCESS\n" : "FAILED\n");
		if (success && !status) success = false;
	}


	return EXIT_SUCCESS;
}
