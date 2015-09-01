## Process this file with automake to produce Makefile.in

AM_CFLAGS =\
	 -Wall\
	 -g

lib_LTLIBRARIES = libclinalg.la

libclinalg_la_SOURCES = \
	clinalg_utils.c\
	matrix/c_matrix.c\
	matrix/c_vector.c\
	matrix/c_vector_int.c\
	matrix/c_matrixops.c\
	matrix/c_vectorops.c\
	linalg/c_linalg_lapack_lu.c\
	linalg/c_linalg_lu.c\
	linalg/c_linalg_lapack_cholesky.c\
	linalg/c_linalg_cholesky.c\
	linalg/c_linalg_lapack_qr.c\
	linalg/c_linalg_qr.c\
	linalg/c_linalg_lapack_sv.c\
	linalg/c_linalg_sv.c\
	linalg/c_linalg_lapack_eigen.c\
	linalg/c_linalg_eigen.c

libclinalg_la_LIBADD = -L$(top_srcdir)/src
libclinalg_la_CPPFLAGS = -I$(top_srcdir)/include\
	-I$(top_srcdir)/include/matrix -I$(top_srcdir)/include/linalg
libclinalg_la_LDFLAGS = -version-info $(CLINALG_LT_VERSION)

