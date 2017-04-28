
#ifndef BLAS_H
#define BLAS_H

#include "settings.h"

#if defined(__GNUC__) && defined(USE_BLAS)
INLINE_PREFIX void
scal(int n, field alpha, field *x, int incx) __attribute__((unused));
INLINE_PREFIX void
axpy(int n, field alpha, const field *x, int incx, field *y, int incy) __attribute__((unused));
INLINE_PREFIX field
dot(int n, const field *x, int incx, const field *y, int incy) __attribute__((unused));
INLINE_PREFIX field
nrm2(int n, const field *x, int incx) __attribute__((unused));
INLINE_PREFIX void
gemv(bool trans, int m, int n, field alpha, const field *A, int ldA,
     const field *x, int incx, field beta, field *y, int incy) __attribute__((unused));
INLINE_PREFIX void
geru(int m, int n, field alpha, const field *x, int incx,
     const field *y, int incy, field *A, int ldA) __attribute__((unused));
INLINE_PREFIX void
gemm(bool transA, bool transB, int m, int n, int k, field alpha,
     const field *A, int ldA, const field *B, int ldB,
     field beta, field *C, int ldC) __attribute__((unused));
#endif

/* ------------------------------------------------------------
 * Scaling a vector
 * ------------------------------------------------------------ */

#ifndef USE_BLAS
HEADER_PREFIX void
scal(int n, field alpha, field *x, int incx);
#else
HEADER_PREFIX void
sscal_(const int *n, const float *alpha, float *x, const int *incx);

HEADER_PREFIX void
dscal_(const int *n, const double *alpha, double *x, const int *incx);

#ifdef USE_COMPLEX
HEADER_PREFIX void
cscal_(const int *n, const float complex *alpha, float complex *x, const int *incx);

HEADER_PREFIX void
zscal_(const int *n, const double complex *alpha, double complex *x, const int *incx);
#endif

INLINE_PREFIX void
scal(int n, field alpha, field *x, int incx)
{
#ifdef USE_FLOAT
#ifdef USE_COMPLEX
  cscal_(&n, &alpha, x, &incx);
#else
  sscal_(&n, &alpha, x, &incx);
#endif
#else
#ifdef USE_COMPLEX
  zscal_(&n, &alpha, x, &incx);
#else
  dscal_(&n, &alpha, x, &incx);
#endif
#endif
}
#endif

/* ------------------------------------------------------------
 * Adding a vector to another vector
 * ------------------------------------------------------------ */

#ifndef USE_BLAS
HEADER_PREFIX void
axpy(int n, field alpha, const field *x, int incx, field *y, int incy);
#else
HEADER_PREFIX void
saxpy_(const int *n, const float *alpha, const float *x, const int *incx, float *y, const int *incy);

HEADER_PREFIX void
daxpy_(const int *n, const double *alpha, const double *x, const int *incx, double *y, const int *incy);

#ifdef USE_COMPLEX
HEADER_PREFIX void
caxpy_(const int *n, const float complex *alpha, const float complex *x, const int *incx, float complex *y, const int *incy);

HEADER_PREFIX void
zaxpy_(const int *n, const double complex *alpha, const double complex *x, const int *incx, double complex *y, const int *incy);
#endif

INLINE_PREFIX void
axpy(int n, field alpha, const field *x, int incx, field *y, int incy)
{
#ifdef USE_FLOAT
#ifdef USE_COMPLEX
  caxpy_(&n, &alpha, x, &incx, y, &incy);
#else
  saxpy_(&n, &alpha, x, &incx, y, &incy);
#endif
#else
#ifdef USE_COMPLEX
  zaxpy_(&n, &alpha, x, &incx, y, &incy);
#else
  daxpy_(&n, &alpha, x, &incy, y, &incy);
#endif
#endif
}
#endif

/* ------------------------------------------------------------
 * Inner product
 * ------------------------------------------------------------ */

#ifndef USE_BLAS
HEADER_PREFIX field
dot(int n, const field *x, int incx, const field *y, int incy);
#else
HEADER_PREFIX float
sdot_(const int *n, const float *x, const int *incx, const float *y, const int *incy);

HEADER_PREFIX double
ddot_(const int *n, const double *x, const int *incx, const double *y, const int *incy);

#ifdef USE_COMPLEX
HEADER_PREFIX float complex
cdot_(const int *n, const float complex *x, const int *incx, const float complex *y, const int *incy);

HEADER_PREFIX double complex
zdot_(const int *n, const double complex *x, const int *incx, const double complex *y, const int *incy);
#endif

INLINE_PREFIX field
dot(int n, const field *x, int incx, const field *y, int incy)
{
#ifdef USE_FLOAT
#ifdef USE_COMPLEX
  return cdot_(&n, x, &incx, y, &incy);
#else
  return sdot_(&n, x, &incx, y, &incy);
#endif
#else
#ifdef USE_COMPLEX
  return zdot_(&n, x, &incx, y, &incy);
#else
  return ddot_(&n, x, &incx, y, &incy);
#endif
#endif
}
#endif

/* ------------------------------------------------------------
 * Euclidean norm
 * ------------------------------------------------------------ */

#ifndef USE_BLAS
HEADER_PREFIX real
nrm2(int n, const field *x, int incx);
#else
HEADER_PREFIX float
snrm2_(const int *n, const float *x, const int *incx);

HEADER_PREFIX double
dnrm2_(const int *n, const double *x, const int *incx);

#ifdef USE_COMPLEX
HEADER_PREFIX float
cnrm2_(const int *n, const float complex *x, const int *incx);

HEADER_PREFIX double
znrm2_(const int *n, const double complex *x, const int *incx);
#endif

INLINE_PREFIX field
nrm2(int n, const field *x, int incx)
{
#ifdef USE_FLOAT
#ifdef USE_COMPLEX
  return cnrm2_(&n, x, &incx);
#else
  return snrm2_(&n, x, &incx);
#endif
#else
#ifdef USE_COMPLEX
  return znrm2_(&n, x, &incx);
#else
  return dnrm2_(&n, x, &incx);
#endif
#endif
}
#endif

/* ------------------------------------------------------------
 * Matrix-vector product
 * ------------------------------------------------------------ */

#ifndef USE_BLAS
HEADER_PREFIX void
gemv(bool trans, int m, int n, field alpha, const field *A, int ldA,
     const field *x, int incx, field beta, field *y, int incy);
#else
HEADER_PREFIX void
sgemv_(const char *trans, const int *m, const int *n,
       const float *alpha, const float *A, const int *ldA,
       const float *x, const int *incx,
       const float *beta, float *y, const int *incy);

HEADER_PREFIX void
dgemv_(const char *trans, const int *m, const int *n,
       const double *alpha, const double *A, const int *ldA,
       const double *x, const int *incx,
       const double *beta, double *y, const int *incy);

#ifdef USE_COMPLEX
HEADER_PREFIX void
cgemv_(const char *trans, const int *m, const int *n,
       const float complex *alpha, const float complex *A, const int *ldA,
       const float complex *x, const int *incx,
       const float complex *beta, float complex *y, const int *incy);

HEADER_PREFIX void
zgemv_(const char *trans, const int *m, const int *n,
       const double complex *alpha, const double complex *A, const int *ldA,
       const double complex *x, const int *incx,
       const double complex *beta, double complex *y, const int *incy);
#endif

INLINE_PREFIX void
gemv(bool trans, int m, int n, field alpha, const field *A, int ldA,
     const field *x, int incx, field beta, field *y, int incy)
{
#ifdef USE_FLOAT
#ifdef USE_COMPLEX
  cgemv_((trans ? "c" : "n"), &m, &n, &alpha, A, &ldA,
	   x, &incx, &beta, y, &incy);
#else
  sgemv_((trans ? "t" : "n"), &m, &n, &alpha, A, &ldA,
	 x, &incx, &beta, y, &incy);
#endif
#else
#ifdef USE_COMPLEX
  zgemv_((trans ? "c" : "n"), &m, &n, &alpha, A, &ldA,
	 x, &incx, &beta, y, &incy);
#else
  dgemv_((trans ? "t" : "n"), &m, &n, &alpha, A, &ldA,
	 x, &incx, &beta, y, &incy);
#endif
#endif
}
#endif

/* ------------------------------------------------------------
 * Rank-one update
 * ------------------------------------------------------------ */

#ifndef USE_BLAS
HEADER_PREFIX void
geru(int m, int n, field alpha, const field *x, int incx,
     const field *y, int incy, field *A, int ldA);
#else
HEADER_PREFIX void
sger_(const int *m, const int *n, const float *alpha,
      const float *x, const int *incx, const float *y, const int *incy,
      float *A, const int *ldA);

HEADER_PREFIX void
dger_(const int *m, const int *n, const double *alpha,
      const double *x, const int *incx, const double *y, const int *incy,
      double *A, const int *ldA);

#ifdef USE_COMPLEX
HEADER_PREFIX void
cgeru_(const int *m, const int *n, const float complex *alpha,
       const float complex *x, const int *incx, const float complex *y, const int *incy,
       float complex *A, const int *ldA);

HEADER_PREFIX void
zgeru_(const int *m, const int *n, const double complex *alpha,
       const double complex *x, const int *incx, const double complex *y, const int *incy,
       double complex *A, const int *ldA);
#endif

INLINE_PREFIX void
geru(int m, int n, field alpha, const field *x, int incx,
     const field *y, int incy, field *A, int ldA)
{
#ifdef USE_FLOAT
#ifdef USE_COMPLEX
  cgeru_(&m, &n, &alpha, x, &incx, y, &incy, A, &ldA);
#else
  sger_(&m, &n, &alpha, x, &incx, y, &incy, A, &ldA);
#endif
#else
#ifdef USE_COMPLEX
  zgeru_(&m, &n, &alpha, x, &incx, y, &incy, A, &ldA);
#else
  dger_(&m, &n, &alpha, x, &incx, y, &incy, A, &ldA);
#endif
#endif
}
#endif

/* ------------------------------------------------------------
 * Matrix multiplication
 * ------------------------------------------------------------ */

#ifndef USE_BLAS
HEADER_PREFIX void
gemm(bool transA, bool transB, int m, int n, int k, field alpha,
     const field *A, int ldA, const field *B, int ldB,
     field beta, field *C, int ldC);
#else
HEADER_PREFIX void
sgemm_(const char *transA, const char *transB,
       const int *m, const int *n, const int *k,
       const float *alpha, const float *A, const int *ldA,
       const float *B, const int *ldB,
       const float *beta, float *C, const int *ldC);

HEADER_PREFIX void
dgemm_(const char *transA, const char *transB,
       const int *m, const int *n, const int *k,
       const double *alpha, const double *A, const int *ldA,
       const double *B, const int *ldB,
       const double *beta, double *C, const int *ldC);

#ifdef USE_COMPLEX
HEADER_PREFIX void
cgemm_(const char *transA, const char *transB,
       const int *m, const int *n, const int *k,
       const float complex *alpha, const float complex *A, const int *ldA,
       const float complex *B, const int *ldB,
       const float complex *beta, float complex *C, const int *ldC);

HEADER_PREFIX void
zgemm_(const char *transA, const char *transB,
       const int *m, const int *n, const int *k,
       const double complex *alpha, const double complex *A, const int *ldA,
       const double complex *B, const int *ldB,
       const double complex *beta, double complex *C, const int *ldC);
#endif

INLINE_PREFIX void
gemm(bool transA, bool transB, int m, int n, int k, field alpha,
     const field *A, int ldA, const field *B, int ldB,
     field beta, field *C, int ldC)
{
#ifdef USE_FLOAT
#ifdef USE_COMPLEX
  cgemm_((transA ? "c" : "n"), (transB ? "c" : "n"),
	 &m, &n, &k, &alpha, A, &ldA, B, &ldB, &beta, C, &ldC);
#else
  sgemm_((transA ? "t" : "n"), (transB ? "t" : "n"),
	 &m, &n, &k, &alpha, A, &ldA, B, &ldB, &beta, C, &ldC);
#endif
#else
#ifdef USE_COMPLEX
  zgemm_((transA ? "c" : "n"), (transB ? "c" : "n"),
	 &m, &n, &k, &alpha, A, &ldA, B, &ldB, &beta, C, &ldC);
#else
  dgemm_((transA ? "t" : "n"), (transB ? "t" : "n"),
	 &m, &n, &k, &alpha, A, &ldA, B, &ldB, &beta, C, &ldC);
#endif
#endif
}
#endif


#endif
