/*
 * private.h
 *
 *  Created on: 2014/10/23
 *      Author: utsugi
 */

#ifndef PRIVATE_H_
#define PRIVATE_H_

/* blas */
#ifdef HAVE_BLAS_H
#include <blas.h>
#else
// level1
extern int		idamax_ (int *n, double *x, int *incx);
extern double	dasum_ (int *n, double *x, int *incx);
extern double	dnrm2_ (int *n, double *x, int *incx);
extern double	ddot_ (int *n, double *x, int *incx, double *y, int *incy);
extern void	dscal_ (int *n, double *alpha, double *x, int *incx);
extern void	daxpy_ (int *n, double *alpha, double *x, int *incx, double *y, int *incy);
extern void	dcopy_ (int *n, double *x, int *incx, double *y, int *incy);
// level2
extern void	dgemv_ (char *trans, int *n, int *m, double *alpha, double *a, int *lda, double *x, int *incx, double *beta, double *y, int *incy);
extern void	dsymv_ (char *uplo, int *n, double *alpha, double *a, int *lda, double *x, int *incx, double *beta, double *y, int *incy);
// level3
extern void	dgemm_ (char *transA, char *transB, int *m, int *n, int *k, double *alpha, double *a, int *lda, double *b, int *ldb, double *beta, double *c, int *ldc);
extern void	dsymm_ (char *side, char *uplo, int *m, int *n, double *alpha, double *a, int *lda, double *b, int *ldb, double *beta, double *c, int *ldc);
#endif

#ifdef HAVE_LAPACK_H
#include <lapack.h>
#else
extern double	dlange_ (char *norm, int *m, int *n, double *data, int *lda, double *w);
extern void	dswap_ (int *n, double *x, int *incx, double *y, int *incy);
extern void	dtrtrs_ (char *uplo, char *trans, char *diag, int *n, int *nrhs, double *a, int *lda, double *b, int *ldb, int *info);

/* LU decomposition */
extern void	dgetrf_ (int *m, int *n, double *a, int *lda, int *ipiv, int *info);
extern void	dgetrs_ (char *trans, int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info);
extern void	dgesv_  (int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info);
extern void	dgetri_ (int *n, double *a, int *lda, int *ipiv, double *work, int *lwork, int *info);

/* Cholesky decomposition */
extern void	dpotrf_ (char *uplo, int *n, double *a, int *lda, int *info);
extern void	dpotrs_ (char *uplo, int *n, int *nrhs, double *a, int *lda, double *b, int *ldb, int *info);
extern void	dpotri_ (char *uplo, int *n, double *a, int *lda, int *info);

/* QR decomposition */
extern void	dgeqrf_ (int *m, int *n, double *data, int *lda, double *tau, double *work, int *lwork, int *info);
extern void	dgeqp3_ (int *m, int *n, double *a, int *lda, int *jpvt, double *tau, double *work, int *lwork, int *info);
extern void	dorgqr_ (int *m, int *n, int *k, double *data, int *lda, double *tau, double *work, int *lwork, int *info);
extern void	dgels_  (char *trans, int *m, int *n, int *nrhs, double *a_data, int *lda, double *b_data, int *ldb, double *w, int *lwork, int *info);
extern void	dgelsy_ (int *m, int *n, int *nrhs, double *a, int *lda, double *b, int *ldb, int *jpvt, double *rcond, int *rank, double *work, int *lwork, int *info);

/* SV decomposition */
extern int		ilaenv_ (int *ispec, char *name, char *opts, int *n1, int *n2, int *n3, int *n4);
extern void	dgesvd_ (char *jobu, char *jobvt, int *m, int *n, double *a_data, int *lda, double *s_data,
						  double *u_data, int *ldu, double *vt_data, int *ldvt, double *w, int *lwork, int *info);
extern void	dgesdd_ (char *jobz, int *m, int *n, double *a_data, int *lda, double *s_data, double *u_data, int *ldu,
						  double *vt_data, int *ldvt, double *work, int *lwork, int *iwork, int *info);
extern void	dgelss_ (int *m, int *n, int *nrhs, double *a_data, int *lda, double *b_data, int *ldb,
						  double *s_data, double *rcond, int *lrank, double *w, int *lwork, int *info);
extern void	dgelsd_ (int *m, int *n, int *nrhs, double *a_data, int *lda, double *b_data, int *ldb,
						  double *s_data, double *rcond, int *lrank, double *w, int *lwork, int *iwork, int *info);
#endif

/* qrupdate: cholinsert/delete */
#ifdef HAVE_QRUPDATE_H
#include <qrupdate.h>
#else
extern void	dch1up_ (int *n, double *L, int *ldr, double *u, double *w);
extern void	dch1dn_ (int *n, double *L, int *ldr, double *u, double *w, int *info);
extern void	dchinx_ (int *n, double *L, int *ldr, int *j, double *u, double *w, int *info);
extern void	dchdex_ (int *n, double *L, int *ldr, int *j, double *w);
extern void	dlup1up_ (int *m, int *n, double *L, int *ldl, double *R, int *ldr, int *p, double *u, double *v, double *w);
extern void	dqr1up_ (int *m, int *n, int *k, double *Q, int *ldq, double *R, int *ldr, double *u, double *v, double *w);
extern void	dqrinc_ (int *m, int *n, int *k, double *Q, int *ldq, double *R, int *ldr, int *j, double *x, double *w);
extern void	dqrinr_ (int *m, int *n, double *Q, int *ldq, double *R, int *ldr, int *j, double *x, double *w);
extern void	dqrdec_ (int *m, int *n, int *k, double *Q, int *ldq, double *R, int *ldr, int *j, double *w);
extern void	dqrder_ (int *m, int *n, double *Q, int *ldq, double *R, int *ldr, int *j, double *w);
#endif


#endif /* PRIVATE_H_ */
