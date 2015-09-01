/*
 * private.h
 *
 *  Created on: 2014/10/23
 *      Author: utsugi
 */

#ifndef PRIVATE_H_
#define PRIVATE_H_

extern const int		izero;
extern const int		ione;
extern const double		dzero;
extern const double		done;
extern const double		dmone;

#ifndef F77CALL
#define F77CALL(CallName) CallName##_
#endif

/* positive / negative infinity  */
#ifndef POSINF
#define POSINF	((+1.)/(+0.))
#endif

#ifndef NEGINF
#define NEGINF	((-1.)/(+0.))
#endif

/* DBL)EPSILOM */
#ifndef DBL_EPSILON
#define DBL_EPSILON		2.2204460492503131e-16
#endif

/* SQRT)DBL)EPSILOM */
#ifndef SQRT_DBL_EPSILON
#define SQRT_DBL_EPSILON	1.4901161193847656e-08
#endif

/* blas */
#ifdef HAVE_BLAS_H
#include <blas.h>
#else
// level1
extern int		F77CALL (idamax) (const int *n, const double *x, const int *incx);
extern double	F77CALL (dasum) (const int *n, const double *x, const int *incx);
extern double	F77CALL (dnrm2) (const int *n, const double *x, const int *incx);
extern double	F77CALL (ddot) (const int *n, const double *x, const int *incx, const double *y, const int *incy);
extern void		F77CALL (dscal) (const int *n, const double *alpha, double *x, const int *incx);
extern void		F77CALL (daxpy) (const int *n, const double *alpha, const double *x, const int *incx, double *y, const int *incy);
extern void		F77CALL (dcopy) (const int *n, const double *x, const int *incx, double *y, const int *incy);
// level2
extern void		F77CALL (dgemv) (const char *trans, const int *n, const int *m,const  double *alpha, const double *a, const int *lda,
						 const double *x, const int *incx, const double *beta, double *y, const int *incy);
extern void		F77CALL (dsymv) (const char *uplo, const int *n, const double *alpha, const double *a, const int *lda, const double *x, const int *incx,
						 const double *beta, double *y, const int *incy);
// level3
extern void		F77CALL (dgemm) (const char *transA, const char *transB, const int *m, const int *n, const int *k, const double *alpha,
								 const double *a, const int *lda, const double *b, const int *ldb, const double *beta, double *c, const int *ldc);
extern void		F77CALL (dsymm) (const char *side, const char *uplo, const int *m, const int *n, const double *alpha, const double *a, const int *lda,
								 const double *b, const int *ldb, const double *beta, double *c, const int *ldc);
#endif

#ifdef HAVE_LAPACK_H
#include <lapack.h>
#else
extern double	F77CALL (dlange) (const char *norm, const int *m, const int *n, const double *data, const int *lda, const double *w);
extern void		F77CALL (dswap) (const int *n, double *x, const int *incx, double *y, const int *incy);
extern void		F77CALL (dtrtrs) (const char *uplo, const char *trans, const char *diag, const int *n, const int *nrhs,
						  	  	  double *a, const int *lda, double *b, const int *ldb, int *info);

/* LU decomposition */
extern void		F77CALL (dgetrf) (const int *m, const int *n, double *a, const int *lda, int *ipiv, int *info);
extern void		F77CALL (dgetrs) (const char *trans, const int *n, const int *nrhs, double *a, const int *lda, int *ipiv, double *b, const int *ldb, int *info);
extern void		F77CALL (dgesv)  (const int *n, const int *nrhs, double *a, const int *lda, int *ipiv, double *b, const int *ldb, int *info);
extern void		F77CALL (dgetri) (const int *n, double *a, const int *lda, int *ipiv, double *work, int *lwork, int *info);

/* Cholesky decomposition */
extern void		F77CALL (dpotrf) (const char *uplo, const int *n, double *a, const int *lda, int *info);
extern void		F77CALL (dpotrs) (const char *uplo, const int *n, const int *nrhs, double *a, const int *lda, double *b, const int *ldb, int *info);
extern void		F77CALL (dpotri) (const char *uplo, const int *n, double *a, const int *lda, int *info);

/* QR decomposition */
extern void		F77CALL (dgeqrf) (const int *m, const int *n, double *data, const int *lda, double *tau, double *work, int *lwork, int *info);
extern void		F77CALL (dgeqp3) (const int *m, const int *n, double *a, const int *lda, int *jpvt, double *tau, double *work, int *lwork, int *info);
extern void		F77CALL (dorgqr) (const int *m, const int *n, const int *k, double *data, const int *lda, double *tau, double *work, int *lwork, int *info);
extern void		F77CALL (dgels)  (const char *trans, const int *m, const int *n, const int *nrhs, double *a, int *lda, double *b, const int *ldb, double *w, int *lwork, int *info);
extern void		F77CALL (dgelsy) (const int *m, const int *n, const int *nrhs, double *a, int *lda, double *b, const int *ldb,
						  	  	  int *jpvt, double *rcond, int *rank, double *work, int *lwork, int *info);

/* SV decomposition */
extern int		F77CALL (ilaenv) (const int *ispec, const char *name, const char *opts, const int *n1, const int *n2, const int *n3, const int *n4);
extern void		F77CALL (dgesvd) (const char *jobu, const char *jobvt, const int *m, const int *n, double *a, const int *lda, double *s,
						  	  	  double *u, const int *ldu, double *vt, const int *ldvt, double *w, int *lwork, int *info);
extern void		F77CALL (dgesdd) (const char *jobz, const int *m, const int *n, double *a, const int *lda, double *s, double *u, const int *ldu,
								  double *vt, const int *ldvt, double *work, int *lwork, int *iwork, int *info);
extern void		F77CALL (dgelss) (const int *m, const int *n, const int *nrhs, double *a, const int *lda, double *b, const int *ldb,
						  	  	  double *s, double *rcond, int *lrank, double *w, int *lwork, int *info);
extern void		F77CALL (dgelsd) (const int *m, const int *n, const int *nrhs, double *a, const int *lda, double *b, const int *ldb,
						  	  	  double *s, double *rcond, int *lrank, double *w, int *lwork, int *iwork, int *info);

/* eigen value decomposition */
extern void		F77CALL (dsyev) (const char *jobz, const char *uplo, const int *n, double *a, const int *lda,
								 double *w, double *work, int *lwork, int *info);
extern void		F77CALL (dsyevd) (const char *jobz, const char *uplo, const int *n, double *a, const int *lda,
								  double *w, double *work, int *lwork, int *iwork, int *liwork, int *info);
extern void		F77CALL (dgeev) (const char *jobvl, const char *jobvr, const int *n, double *a, const int *lda,
								 double *wr, double *wi, double *vl, const int *ldvl, double *vr, const int *ldvr,
								 double *work, int *lwork, int *info);
extern void		F77CALL (dgees) (const char *jobvs, const char *sort, void *select, const int *n, double *a, const int *lda, int *sdim,
				  	  	  	  	 double *wr, double *wi, double *vs, const int *ldvs, double *work, int *lwork, int *bwork, int *info);

#endif

/* qrupdate: cholinsert/delete */
#ifdef HAVE_QRUPDATE_H
#include <qrupdate.h>
#else
extern void		F77CALL (dch1up) (const int *n, double *L, const int *ldr, double *u, double *w);
extern void		F77CALL (dch1dn) (const int *n, double *L, const int *ldr, double *u, double *w, int *info);
extern void		F77CALL (dchinx) (const int *n, double *L, const int *ldr, const int *j, double *u, double *w, int *info);
extern void		F77CALL (dchdex) (const int *n, double *L, const int *ldr, const int *j, double *w);
extern void		F77CALL (dlup1up) (const int *m, const int *n, double *L, const int *ldl, double *R, const int *ldr, int *p, double *u, double *v, double *w);
extern void		F77CALL (dqr1up) (const int *m, const int *n, const int *k, double *Q, const int *ldq, double *R, const int *ldr, double *u, double *v, double *w);
extern void		F77CALL (dqrinc) (const int *m, const int *n, const int *k, double *Q, const int *ldq, double *R, const int *ldr, const int *j, double *x, double *w);
extern void		F77CALL (dqrinr) (const int *m, const int *n, double *Q, const int *ldq, double *R, const int *ldr, const int *j, double *x, double *w);
extern void		F77CALL (dqrdec) (const int *m, const int *n, const int *k, double *Q, const int *ldq, double *R, const int *ldr, const int *j, double *w);
extern void		F77CALL (dqrder) (const int *m, const int *n, double *Q, const int *ldq, double *R, const int *ldr, const int *j, double *w);
#endif


#endif /* PRIVATE_H_ */
