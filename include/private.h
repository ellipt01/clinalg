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
extern int		idamax_ (const int *n, const double *x, const int *incx);
extern double	dasum_ (const int *n, const double *x, const int *incx);
extern double	dnrm2_ (const int *n, const double *x, const int *incx);
extern double	ddot_ (const int *n, const double *x, const int *incx, const double *y, const int *incy);
extern void	dscal_ (const int *n, const double *alpha, double *x, const int *incx);
extern void	daxpy_ (const int *n, const double *alpha, const double *x, const int *incx, double *y, const int *incy);
extern void	dcopy_ (const int *n, const double *x, const int *incx, double *y, const int *incy);
// level2
extern void	dgemv_ (const char *trans, const int *n, const int *m,const  double *alpha, const double *a, const int *lda,
						 const double *x, const int *incx, const double *beta, double *y, const int *incy);
extern void	dsymv_ (const char *uplo, const int *n, const double *alpha, const double *a, const int *lda, const double *x, const int *incx,
						 const double *beta, double *y, const int *incy);
// level3
extern void	dgemm_ (const char *transA, const char *transB, const int *m, const int *n, const int *k, const double *alpha,
						 const double *a, const int *lda, const double *b, const int *ldb, const double *beta, double *c, const int *ldc);
extern void	dsymm_ (const char *side, const char *uplo, const int *m, const int *n, const double *alpha, const double *a, const int *lda,
						 const double *b, const int *ldb, const double *beta, double *c, const int *ldc);
#endif

#ifdef HAVE_LAPACK_H
#include <lapack.h>
#else
extern double	dlange_ (const char *norm, const int *m, const int *n, const double *data, const int *lda, const double *w);
extern void	dswap_ (const int *n, double *x, const int *incx, double *y, const int *incy);
extern void	dtrtrs_ (const char *uplo, const char *trans, const char *diag, const int *n, const int *nrhs,
						  double *a, const int *lda, double *b, const int *ldb, int *info);

/* LU decomposition */
extern void	dgetrf_ (const int *m, const int *n, double *a, const int *lda, int *ipiv, int *info);
extern void	dgetrs_ (const char *trans, const int *n, const int *nrhs, double *a, const int *lda, int *ipiv, double *b, const int *ldb, int *info);
extern void	dgesv_  (const int *n, const int *nrhs, double *a, const int *lda, int *ipiv, double *b, const int *ldb, int *info);
extern void	dgetri_ (const int *n, double *a, const int *lda, int *ipiv, double *work, int *lwork, int *info);

/* Cholesky decomposition */
extern void	dpotrf_ (const char *uplo, const int *n, double *a, const int *lda, int *info);
extern void	dpotrs_ (const char *uplo, const int *n, const int *nrhs, double *a, const int *lda, double *b, const int *ldb, int *info);
extern void	dpotri_ (const char *uplo, const int *n, double *a, const int *lda, int *info);

/* QR decomposition */
extern void	dgeqrf_ (const int *m, const int *n, double *data, const int *lda, double *tau, double *work, int *lwork, int *info);
extern void	dgeqp3_ (const int *m, const int *n, double *a, const int *lda, int *jpvt, double *tau, double *work, int *lwork, int *info);
extern void	dorgqr_ (const int *m, const int *n, const int *k, double *data, const int *lda, double *tau, double *work, int *lwork, int *info);
extern void	dgels_  (const char *trans, const int *m, const int *n, const int *nrhs, double *a_data, int *lda, double *b_data, const int *ldb, double *w, int *lwork, int *info);
extern void	dgelsy_ (const int *m, const int *n, const int *nrhs, double *a, int *lda, double *b, const int *ldb,
						  int *jpvt, double *rcond, int *rank, double *work, int *lwork, int *info);

/* SV decomposition */
extern int		ilaenv_ (const int *ispec, const char *name, const char *opts, const int *n1, const int *n2, const int *n3, const int *n4);
extern void	dgesvd_ (const char *jobu, const char *jobvt, const int *m, const int *n, double *a_data, const int *lda, double *s_data,
						  double *u_data, const int *ldu, double *vt_data, const int *ldvt, double *w, int *lwork, int *info);
extern void	dgesdd_ (const char *jobz, const int *m, const int *n, double *a_data, const int *lda, double *s_data, double *u_data, const int *ldu,
						  double *vt_data, const int *ldvt, double *work, int *lwork, int *iwork, int *info);
extern void	dgelss_ (const int *m, const int *n, const int *nrhs, double *a_data, const int *lda, double *b_data, const int *ldb,
						  double *s_data, double *rcond, int *lrank, double *w, int *lwork, int *info);
extern void	dgelsd_ (const int *m, const int *n, const int *nrhs, double *a_data, const int *lda, double *b_data, const int *ldb,
						  double *s_data, double *rcond, int *lrank, double *w, int *lwork, int *iwork, int *info);
#endif

/* qrupdate: cholinsert/delete */
#ifdef HAVE_QRUPDATE_H
#include <qrupdate.h>
#else
extern void	dch1up_ (const int *n, double *L, const int *ldr, double *u, double *w);
extern void	dch1dn_ (const int *n, double *L, const int *ldr, double *u, double *w, int *info);
extern void	dchinx_ (const int *n, double *L, const int *ldr, const int *j, double *u, double *w, int *info);
extern void	dchdex_ (const int *n, double *L, const int *ldr, const int *j, double *w);
extern void	dlup1up_ (const int *m, const int *n, double *L, const int *ldl, double *R, const int *ldr, int *p, double *u, double *v, double *w);
extern void	dqr1up_ (const int *m, const int *n, const int *k, double *Q, const int *ldq, double *R, const int *ldr, double *u, double *v, double *w);
extern void	dqrinc_ (const int *m, const int *n, const int *k, double *Q, const int *ldq, double *R, const int *ldr, const int *j, double *x, double *w);
extern void	dqrinr_ (const int *m, const int *n, double *Q, const int *ldq, double *R, const int *ldr, const int *j, double *x, double *w);
extern void	dqrdec_ (const int *m, const int *n, const int *k, double *Q, const int *ldq, double *R, const int *ldr, const int *j, double *w);
extern void	dqrder_ (const int *m, const int *n, double *Q, const int *ldq, double *R, const int *ldr, const int *j, double *w);
#endif


#endif /* PRIVATE_H_ */
