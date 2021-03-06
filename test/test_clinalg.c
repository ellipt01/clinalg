/*
 * test.c
 *
 *  Created on: 2014/04/14
 *      Author: utsugi
 */

#include <stdlib.h>
#include <clinalg.h>
#include "test_clinalg.h"

int		size1;
int		size2;

const int	n_test_matrix_func = 2;

const char	*test_matrix_func_name[] = {
	"test_matrix_mv",
	"test_matrix_mm"
};

bool	(*test_matrix_func_ptr[]) (void) = {
		test_matrix_mv,
		test_matrix_mm
};

const int	n_test_cholesky_func = 6;

const char	*test_cholesky_func_name[] = {
	"test_cholesky_decomp",
	"test_cholesky_svx   ",
	"test_cholesky_1up   ",
	"test_cholesky_1down ",
	"test_cholesky_insert",
	"test_cholesky_delete"
};

bool	(*test_cholesky_func_ptr[]) (void) = {
	test_cholesky_decomp,
	test_cholesky_svx,
	test_cholesky_1up,
	test_cholesky_1down,
	test_cholesky_insert,
	test_cholesky_delete
};

const int	n_test_LU_func = 5;

const char	*test_LU_func_name[] = {
	"test_LU_decomp      ",
	"test_LU_solve       ",
	"test_LU_svx         ",
	"test_LU_invert      ",
	"test_LU_1up         "
};

bool	(*test_LU_func_ptr[]) (void) = {
	test_LU_decomp,
	test_LU_solve,
	test_LU_svx,
	test_LU_invert,
	test_LU_1up
};

const int	n_test_QR_func = 10;

const char	*test_QR_func_name[] = {
	"test_QR_decomp      ",
	"test_QR_decomp_econ ",
	"test_QR_solve       ",
	"test_lsQ_solve      ",
	"test_QR_Rsolve      ",
	"test_QR_1up         ",
	"test_QR_colinsert   ",
	"test_QR_rowinsert   ",
	"test_QR_coldelete   ",
	"test_QR_rowdelete   "
};

bool	(*test_QR_func_ptr[]) (void) = {
	test_QR_decomp,
	test_QR_decomp_econ,
	test_QR_solve,
	test_lsQ_solve,
	test_QR_Rsolve,
	test_QR_1up,
	test_QR_colinsert,
	test_QR_rowinsert,
	test_QR_coldelete,
	test_QR_rowdelete
};

const int	n_test_SV_func = 3;

const char	*test_SV_func_name[] = {
	"test_SV_decomp      ",
	"test_SV_solve       ",
	"test_SV_lsd_solve   "
};

bool	(*test_SV_func_ptr[]) (void) = {
	test_SV_decomp,
	test_SV_solve,
	test_SV_lsd_solve
};

#include <time.h>

int
main (void)
{
	int		i;
	bool	success = true;

	srand (time (NULL));

	size1 = (int) 100 * rand () / RAND_MAX + 300;
	size2 = (int) 100 * rand () / RAND_MAX + 300;
	fprintf (stderr, "size1 = %d, size2 = %d\n\n", (int) size1, (int) size2);

	fprintf (stderr, "*** test_matrix ***\n");
	for (i = 0; i < n_test_matrix_func; i++) {
		bool	status = test_matrix_func_ptr[i] ();
		fprintf (stderr, "\t %s ...", test_matrix_func_name[i]);
		fprintf (stderr, status ? "SUCCESS\n" : "FAILED\n");
		if (success && !status) success = false;
	}

	fprintf (stderr, "\n*** test_cholesky ***\n");
	for (i = 0; i < n_test_cholesky_func; i++) {
		bool	status = test_cholesky_func_ptr[i] ();
		fprintf (stderr, "\t %s ...", test_cholesky_func_name[i]);
		fprintf (stderr, status ? "SUCCESS\n" : "FAILED\n");
		if (success && !status) success = false;
	}

	fprintf (stderr, "\n*** test_LU ***\n");
	for (i = 0; i < n_test_LU_func; i++) {
		bool	status = test_LU_func_ptr[i] ();
		fprintf (stderr, "\t %s ...", test_LU_func_name[i]);
		fprintf (stderr, status ? "SUCCESS\n" : "FAILED\n");
		if (success && !status) success = false;
	}

	fprintf (stderr, "\n*** test_QR ***\n");
	for (i = 0; i < n_test_QR_func; i++) {
		bool	status = test_QR_func_ptr[i] ();
		fprintf (stderr, "\t %s ...", test_QR_func_name[i]);
		fprintf (stderr, status ? "SUCCESS\n" : "FAILED\n");
		if (success && !status) success = false;
	}

	fprintf (stderr, "\n*** test_SV ***\n");
	for (i = 0; i < n_test_SV_func; i++) {
		bool	status = test_SV_func_ptr[i] ();
		fprintf (stderr, "\t %s ...", test_SV_func_name[i]);
		fprintf (stderr, status ? "SUCCESS\n" : "FAILED\n");
		if (success && !status) success = false;
	}

	return EXIT_SUCCESS;
}
