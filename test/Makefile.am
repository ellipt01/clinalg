#######################################
# The list of executables we are building seperated by spaces
# the 'bin_' indicates that these build products will be installed
# in the $(bindir) directory. For test_clinalg /usr/bin
#bin_PROGRAMS=test_clinalgProgram

# Because a.out is only a sample program we don't want it to be installed.
# The 'noinst_' prefix indicates that the following targets are not to be
# installed.
noinst_PROGRAMS=test_clinalg

#######################################
# Build information for each executable. The variable name is derived
# by use the name of the executable with each non alpha-numeric character is
# replaced by '_'. So a.out becomes a_out and the appropriate suffex added.
# '_SOURCES' for test_clinalg.

# Sources for the a.out 
test_clinalg_SOURCES=\
	test_cholesky.c\
	test_qr.c\
	test_utils.c\
	test_clinalg.c

# Linker options for a.out
test_clinalg_LDFLAGS = $(top_srcdir)/src/libclinalg.la

# Compiler options for a.out
test_clinalg_CPPFLAGS =\
	-I$(top_srcdir)/include\
	-I$(top_srcdir)/include/matrix\
	-I$(top_srcdir)/include/linalg\
	-I$(top_srcdir)/test 
