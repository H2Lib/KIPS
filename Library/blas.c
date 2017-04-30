
#include "blas.h"

#include "basic.h"

#ifndef USE_BLAS

void
scal(int n, field alpha, field *x, int incx)
{
  int i;

  if(incx == 1)
    for(i=0; i<n; i++)
      x[i] *= alpha;
  else
    for(i=0; i<n; i++)
      x[i*incx] *= alpha;
}

void
axpy(int n, field alpha, const field *x, int incx, field *y, int incy)
{
  int i;

  if(incx == 1 && incy == 1)
    for(i=0; i<n; i++)
      y[i] += alpha * x[i];
  else
    for(i=0; i<n; i++)
      y[i*incy] += alpha * x[i*incx];
}

field
dot(bool conjx, bool conjy, int n,
    const field *x, int incx, const field *y, int incy)
{
  field sum;
  int i;

  sum = 0.0;

  if(incx == 1 && incy == 1) {
    if(conjx) {
      if(conjy)
	for(i=0; i<n; i++)
	  sum += CONJ(x[i]) * CONJ(y[i]);
      else
	for(i=0; i<n; i++)
	  sum += CONJ(x[i]) * y[i];
    }
    else {
      if(conjy)
	for(i=0; i<n; i++)
	  sum += x[i] * CONJ(y[i]);
      else
	for(i=0; i<n; i++)
	  sum += x[i] * y[i];
    }
  }
  else {
    if(conjx) {
      if(conjy)
	for(i=0; i<n; i++)
	  sum += CONJ(x[i*incx]) * CONJ(y[i*incx]);
      else
	for(i=0; i<n; i++)
	  sum += CONJ(x[i*incx]) * y[i*incx];
    }
    else {
      if(conjy)
	for(i=0; i<n; i++)
	  sum += x[i*incx] * CONJ(y[i*incx]);
      else
	for(i=0; i<n; i++)
	  sum += x[i*incx] * y[i*incx];
    }
  }

  return sum;
}

real
nrm2(int n, const field *x, int incx)
{
  real sum;
  int i;

  sum = 0.0;

  if(incx == 1 && incy == 1)
    for(i=0; i<n; i++)
      sum += ABSSQR(x[i]);
  else
    for(i=0; i<n; i++)
      sum += ABSSQR(x[i*incx]);

  return REAL_SQRT(sum);
}

void
gemv(bool trans, int m, int n, field alpha, const field *A, int ldA,
     const field *x, int incx, field beta, field *y, int incy)
{
  int j;

  if(trans) {
    if(beta != 1.0)
      scal(n, beta, y, incy);

    for(j=0; j<n; j++)
      y[j] += dot(true, false, m, A+j*ldA, 1, x, incx);
  }
  else {
    if(beta != 1.0)
      scal(m, beta, y, incy);

    for(j=0; j<n; j++)
      axpy(m, x[j], A+j*ldA, 1, y, incy);
  }
}

void
ger(bool conjy, int m, int n, field alpha, const field *x, int incx,
    const field *y, int incy, field *A, int ldA)
{
  int j;

  if(conjy) {
    for(j=0; j<n; j++)
      axpy(m, CONJ(y[j]), x, 1, A+j*ldA, 1);
  }
  else {
    for(j=0; j<n; j++)
      axpy(m, y[j], x, 1, A+j*ldA, 1);
  }
}

void
gemm(bool transA, bool transB, int m, int n, int k, field alpha,
     const field *A, int ldA, const field *B, int ldB,
     field beta, field *C, int ldC)
{
  field sum;
  int i, j;

  if(transB) {
    if(transA) {
      for(j=0; j<n; j++) {
	if(beta != 1.0)
	  scal(m, beta, C+j*ldC, 1);

	for(i=0; i<m; i++)
	  C[i+j*ldC] += alpha * dot(true, true, k, A+i*ldA, 1, B+j, ldB);
      }
    }
    else {
      if(beta != 1.0)
	for(j=0; j<n; j++)
	  scal(m, beta, C+j*ldC, 1);

      for(j=0; j<k; j++)
	ger(true, m, n, alpha, A+j*ldA, 1, B+j*ldB, 1, C, ldC);
    }
  }
  else {
    for(j=0; j<n; j++) {
      if(beta != 1.0)
	scal(m, beta, C+j*ldC, 1);

      gemv(transA, m, k, alpha, A, ldA, B+j*ldB, 1, C+j*ldC, 1);
    }
  }
}

#endif
