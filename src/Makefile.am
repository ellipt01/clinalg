## Process this file with automake to produce Makefile.in

INCLUDES =\
	 -I$(top_srcdir)/include\
	 -I$(top_srcdir)/include/matrix

AM_CFLAGS =\
	 -Wall\
	 -g

lib_LTLIBRARIES = libclinalg.la

libclinalg_la_SOURCES = \
	c_linalg_utils.c\
	matrix/c_matrix.c\
	matrix/c_vector.c\
	matrix/c_vector_int.c\
	matrix/c_matrixops.c\
	matrix/c_vectorops.c\
	c_linalg_lu.c\
	c_linalg_cholesky.c\
	c_linalg_qr.c\
	c_linalg_sv.c

libclinalg_la_CPPFLAGS =\
	-I$(top_srcdir)/include/matrix -I$(top_srcdir)/include/linalg

libclinalg_la_LDFLAGS = 

libclinalg_la_LIBADD = -lm

