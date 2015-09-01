#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])
AC_INIT([clinalg], [0.2.2])
AC_CONFIG_HEADERS([config.h])

CLINALG_CURRENT=2
CLINALG_REVISION=0
CLINALG_AGE=1
CLINALG_LT_VERSION="${CLINALG_CURRENT}:${CLINALG_REVISION}:${CLINALG_AGE}"
AC_SUBST(CLINALG_LT_VERSION)

AM_INIT_AUTOMAKE([subdir-objects -Wall])

m4_ifdef([AM_PROG_AR], [AM_PROG_AR([ACT-IF-FAIL])])

# Checks for programs.
AC_PROG_CC
AM_PROG_CC_C_O
AC_PROG_CPP
AC_PROG_CXX
AC_PROG_AWK
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PROG_RANLIB
AC_PROG_LIBTOOL

# Checks for libraries.
AC_CHECK_LIB([m], [floor], [])
AC_CHECK_LIB([gfortran], [main],[])

AC_CHECK_LIB(
	[qrupdate],
	[dchinx_],
	[],
	[AC_MSG_ERROR("ERROR: cannot find library for -lqrupdate. check LDFLAGS and retry - stop.")]
)
AC_CHECK_FUNCS([dchinx_ dchdex_])

AC_ARG_WITH(
	[atlas],
	[AC_HELP_STRING([--with-atlas[[=DIR]]],
	[ATLAS library directory])],
	[atlas_libdir="$withval"],
	[atlas_libdir=""])
if test "$atlas_libdir"; then
	LDFLAGS="$LDFLAGS -L$atlas_libdir"
	use_atlas=yes
fi

AC_ARG_WITH(
	[tatlas],
	[AC_HELP_STRING([--with-tatlas[[=DIR]]],
	[parallel ATLAS library directory])],
	[tatlas_libdir="$withval"],
	[tatlas_libdir=""])
if test "$tatlas_libdir"; then
	LDFLAGS="$LDFLAGS -L$tatlas_libdir"
	use_tatlas=yes
fi

AC_ARG_WITH(
	[satlas],
	[AC_HELP_STRING([--with-satlas[[=DIR]]],
	[serial ATLAS library directory])],
	[satlas_libdir="$withval"],
	[satlas_libdir=""])
if test "$satlas_libdir"; then
	LDFLAGS="$LDFLAGS -L$satlas_libdir"
	use_satlas=yes
fi

AC_ARG_WITH(
	[atlas-incdir],
	[AC_HELP_STRING([--with-atlas-incdir[[=DIR]]],
	[ATLAS include directory])],
	[atlas_incdir="$withval"],
	[atlas_incdir=""])
if test "$atlas_incdir"; then
	CPPFLAGS="$CPPFLAGS -I$atlas_incdir"
fi

if test "x$use_atlas" == xyes; then
	AC_CHECK_LIB([f77blas],[main],[])
elif test "x$use_tatlas" == xyes; then
	AC_CHECK_LIB([tatlas],[main],[])
elif test "x$use_satlas" == xyes; then
	AC_CHECK_LIB([satlas],[main],[])
else
	AC_CHECK_LIB([blas],[main],[])
fi

# Checks for header files.
AC_CHECK_HEADER([stdlib.h])
AC_CHECK_HEADER([string.h])
AC_CHECK_HEADER([blas.h])
AC_CHECK_HEADER([lapack.h])
AC_CHECK_HEADER([qrupdate.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([floor pow sqrt log strrchr])
AC_CHECK_FUNCS([dpotrf_ dpotrs_])

# pkgconfig
CLINALG_DESCRIPTION='C interface of lapack, blas and qrupdate'
AC_SUBST([CLINALG_DESCRIPTION])
AC_CONFIG_FILES([clinalg.pc:pkgconfig/clinalg.pc.in])

AC_CONFIG_FILES([Makefile
                 include/Makefile
                 src/Makefile
                 test/Makefile])
AC_OUTPUT
