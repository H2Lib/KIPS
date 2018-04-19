
/* ------------------------------------------------------------
 * This is the file "amatrix.c" of the KIPS package.
 * All rights reserved, Steffen Boerm 2009
 * ------------------------------------------------------------ */

#include "amatrix.h"

#include "settings.h"
#include "basic.h"

#include <math.h>
#include <stdio.h>

static uint active_amatrix = 0;

/* ------------------------------------------------------------
 * Constructors and destructors
 * ------------------------------------------------------------ */

pamatrix
init_amatrix(pamatrix A, uint rows, uint cols)
{
  assert(A != NULL);

  A->a = (rows > 0 && cols > 0 ? allocmatrix(rows, cols) : NULL);
  A->ld = rows;
  A->rows = rows;
  A->cols = cols;
  A->owner = NULL;

#ifdef USE_OPENMP
#pragma omp atomic
#endif
  active_amatrix++;

  return A;
}

pamatrix
init_sub_amatrix(pamatrix A, pamatrix src, uint rows, uint roff,
		 uint cols, uint coff)
{
  longindex ld = src->ld;

  assert(A != NULL);
  assert(src != NULL);
  assert(roff + rows <= src->rows);
  assert(coff + cols <= src->cols);

  A->a = src->a + roff + ld * coff;
  A->ld = src->ld;
  A->rows = rows;
  A->cols = cols;
  A->owner = src;

#ifdef USE_OPENMP
#pragma omp atomic
#endif
  active_amatrix++;

  return A;
}

void
uninit_amatrix(pamatrix A)
{
  assert(A != NULL);

  if (!A->owner && A->a != NULL) {
    freemem(A->a);
  }

  assert(active_amatrix > 0);

#ifdef USE_OPENMP
#pragma omp atomic
#endif
  active_amatrix--;
}

pamatrix
new_amatrix(uint rows, uint cols)
{
  pamatrix  A;

  A = (pamatrix) allocmem(sizeof(amatrix));

  init_amatrix(A, rows, cols);

  return A;
}

pamatrix
new_sub_amatrix(pamatrix src, uint rows, uint roff, uint cols, uint coff)
{
  pamatrix  A;

  A = (pamatrix) allocmem(sizeof(amatrix));

  init_sub_amatrix(A, src, rows, roff, cols, coff);

  return A;
}

void
del_amatrix(pamatrix A)
{
  assert(A != NULL);

  uninit_amatrix(A);
  freemem(A);
}

void
resize_amatrix(pamatrix A, uint rows, uint cols)
{
  assert(A != NULL);
  assert(A->owner == NULL);

  if (rows != A->rows || cols != A->cols) {
    if (A->a != NULL) {
      freemem(A->a);
    }
    A->a = allocmatrix(rows, cols);
    A->rows = rows;
    A->cols = cols;
    A->ld = rows;
  }
}

void
resizecopy_amatrix(pamatrix A, uint rows, uint cols)
{
  pfield    new_a;
  longindex ldA, ldN;
  uint      i, j;

  assert(A != NULL);
  assert(A->owner == NULL);

  if (rows != A->rows || cols != A->cols) {
    new_a = allocmatrix(rows, cols);
    ldA = A->ld;
    ldN = rows;

    for (j = 0; j < cols && j < A->cols; j++) {
      for (i = 0; i < rows && i < A->rows; i++)
	new_a[i + j * ldN] = A->a[i + j * ldA];
      for (; i < rows; i++)
	new_a[i + j * ldN] = 0.0;
    }
    for (; j < cols; j++)
      for (i = 0; i < rows; i++)
	new_a[i + j * ldN] = 0.0;

    if (A->a != NULL) {
      freemem(A->a);
    }

    A->a = new_a;
    A->rows = rows;
    A->cols = cols;
    A->ld = rows;
  }
}

/* ------------------------------------------------------------
 * Statistics
 * ------------------------------------------------------------ */

uint
getactives_amatrix()
{
  return active_amatrix;
}

size_t
getsize_amatrix(pcamatrix A)
{
  size_t    sz;

  sz = sizeof(amatrix);
  if (A->owner == NULL)
    sz += (size_t) sizeof(field) * A->rows * A->cols;

  return sz;
}

size_t
getsize_heap_amatrix(pcamatrix A)
{
  size_t    sz;

  sz = 0;
  if (A->owner == NULL)
    sz += (size_t) sizeof(field) * A->rows * A->cols;

  return sz;
}

/* ------------------------------------------------------------
 * Simple utility functions
 * ------------------------------------------------------------ */

void
clear_amatrix(pamatrix A)
{
  uint      rows = A->rows;
  longindex ldA = A->ld;
  uint      i, j;

  for (j = 0; j < A->cols; j++)
    for (i = 0; i < rows; i++)
      A->a[i + j * ldA] = 0.0;
}

void
clear_lower_amatrix(bool strict, pamatrix A)
{
  uint      rows = A->rows;
  uint      cols = A->cols;
  longindex ldA = A->ld;
  uint      i, j;

  if (strict) {
    for (j = 0; j < cols; j++) {
      for (i = j + 1; i < rows; i++) {
	A->a[i + j * ldA] = 0.0;
      }
    }
  }
  else {
    for (j = 0; j < cols; j++) {
      for (i = j; i < rows; i++) {
	A->a[i + j * ldA] = 0.0;
      }
    }
  }
}

void
clear_upper_amatrix(bool strict, pamatrix A)
{
  uint      rows = A->rows;
  uint      cols = A->cols;
  longindex ldA = A->ld;
  uint      i, j;

  if (strict) {
    for (j = 0; j < cols; j++) {
      for (i = 0; i < UINT_MIN(j, rows); i++) {
	A->a[i + j * ldA] = 0.0;
      }
    }
  }
  else {
    for (j = 0; j < cols; j++) {
      for (i = 0; i <= UINT_MIN(j, rows - 1); i++) {
	A->a[i + j * ldA] = 0.0;
      }
    }
  }
}

void
identity_amatrix(pamatrix A)
{
  uint      rows = A->rows;
  uint      cols = A->cols;
  longindex ldA = A->ld;
  uint      i, j;

  for (j = 0; j < cols; j++) {
    for (i = 0; i < rows; i++) {
      A->a[i + j * ldA] = 0.0;
    }
  }

  for (i = 0; i < cols && i < rows; i++) {
    A->a[i + i * ldA] = 1.0;
  }
}

void
random_amatrix(pamatrix A)
{
  uint      rows = A->rows;
  uint      cols = A->cols;
  longindex ldA = A->ld;

  uint      i, j;

  for (j = 0; j < cols; j++) {
    for (i = 0; i < rows; i++) {
      A->a[i + j * ldA] = FIELD_RAND();
    }
  }
}

void
random_invertible_amatrix(real alpha, pamatrix A)
{
  pfield    Aa = A->a;
  uint      rows = A->rows;
  longindex ldA = A->ld;
  real      sum;
  uint      i, j;

  assert(rows == A->cols);

  for (j = 0; j < rows; j++) {
    for (i = 0; i < rows; i++) {
      Aa[i + j * ldA] = FIELD_RAND();
    }
    Aa[j + j * ldA] = 0.0;
  }

  /* Ensure strict diagonal dominance if alpha > 1 */
  for (j = 0; j < rows; j++) {
    sum = 0.0;
    for (i = 0; i < rows; i++)
      sum += ABS(Aa[i + j * ldA]);
    Aa[j + j * ldA] = (sum == 0.0 ? alpha : alpha * sum);
  }
}

void
random_selfadjoint_amatrix(pamatrix A)
{
  pfield    Aa = A->a;
  uint      rows = A->rows;
  longindex ldA = A->ld;
  uint      i, j;

  assert(rows == A->cols);

  for (j = 0; j < rows; j++) {
    for (i = 0; i < j; i++) {
      Aa[i + j * ldA] = FIELD_RAND();
      Aa[j + i * ldA] = CONJ(Aa[i + j * ldA]);
    }
    Aa[j + j * ldA] = REAL_RAND();
  }
}

void
random_spd_amatrix(real alpha, pamatrix A)
{
  pfield    Aa = A->a;
  uint      rows = A->rows;
  longindex ldA = A->ld;
  real      sum;
  uint      i, j;

  assert(rows == A->cols);

  for (j = 0; j < rows; j++) {
    for (i = 0; i < j; i++) {
      Aa[i + j * ldA] = FIELD_RAND();
      Aa[j + i * ldA] = CONJ(Aa[i + j * ldA]);
    }
    Aa[j + j * ldA] = 0.0;
  }

  /* Ensure strict diagonal dominance if alpha > 0 */
  for (j = 0; j < rows; j++) {
    sum = 0.0;
    for (i = 0; i < rows; i++)
      sum += ABS(Aa[i + j * ldA]);
    Aa[j + j * ldA] = (sum == 0.0 ? alpha : alpha * sum);
  }
}

void
copy_amatrix(bool transA, pcamatrix A, pamatrix B)
{
  if (transA) {
    assert(A->rows == B->cols);
    assert(A->cols == B->rows);

    copy_sub_amatrix(true, A, B);
  }
  else {
    assert(A->rows == B->rows);
    assert(A->cols == B->cols);

    copy_sub_amatrix(false, A, B);
  }
}

pamatrix
clone_amatrix(pcamatrix src)
{
  pamatrix  A;

  A = new_amatrix(src->rows, src->cols);
  copy_amatrix(false, src, A);

  return A;
}

void
copy_sub_amatrix(bool transA, pcamatrix A, pamatrix B)
{
  pcfield Aa = A->a;
  pfield Ba = B->a;
  longindex ldA = A->ld;
  longindex ldB = B->ld;
  uint      i, j, rows, cols;

  if (transA) {
    rows = UINT_MIN(A->cols, B->rows);
    cols = UINT_MIN(A->rows, B->cols);

    for (i = 0; i < rows; i++) {
      for (j = 0; j < cols; j++) {
	Ba[i + j * ldB] = CONJ(Aa[j + i * ldA]);
      }
    }
  }
  else {
    rows = UINT_MIN(A->rows, B->rows);
    cols = UINT_MIN(A->cols, B->cols);

    for (j = 0; j < cols; j++) {
      for (i = 0; i < rows; i++) {
	Ba[i + j * ldB] = Aa[i + j * ldA];
      }
    }
  }
}

void
print_amatrix(pcamatrix A)
{
  uint      rows = A->rows;
  uint      cols = A->cols;
  longindex ldA = A->ld;
  uint      i, j;

  (void) printf("amatrix(%u,%u,%u)\n", rows, cols, A->ld);
  if (rows == 0 || cols == 0)
    return;

#ifdef USE_COMPLEX
  for (i = 0; i < rows; i++) {
    (void) printf("  (%f+%fi", REAL(A->a[i]), IMAG(A->a[i]));
    for (j = 1; j < cols; j++)
      (void) printf(" | %f+%fi", REAL(A->a[i+j*ldA]), IMAG(A->a[i+j*ldA]));
    (void) printf(")\n");
  }
#else
  for (i = 0; i < rows; i++) {
    (void) printf("  (%f", A->a[i]);
    for (j = 1; j < cols; j++)
      (void) printf(" | %f", A->a[i+j*ldA]);
    (void) printf(")\n");
  }
#endif
}

void
print_matlab_amatrix(pcamatrix A)
{
  uint      rows = A->rows;
  uint      cols = A->cols;
  longindex ldA = A->ld;
  uint      i, j;

  if (rows == 0)
    (void) printf("  [ ]\n");
  else if (rows == 1) {
    if (cols == 0)
      (void) printf("  [ ]\n");
    else {
#ifdef USE_COMPLEX
      (void) printf("  [%f+%fi", REAL(A->a[0]), IMAG(A->a[0]));
      for (j = 1; j < cols; j++)
	(void) printf(" %f+%fi", REAL(A->a[j*ldA]), IMAG(A->a[j*ldA]));
      (void) printf("]\n");
#else
      (void) printf("  [%f", A->a[0]);
      for (j = 1; j < cols; j++)
	(void) printf(" %f", A->a[j*ldA]);
      (void) printf("]\n");
#endif
    }
  }
  else {
    if (cols == 0) {
      (void) printf("  [");
      for (i = 1; i < rows; i++)
	(void) printf(" ;");
      (void) printf(" ]\n");
    }
    else {
#ifdef USE_COMPLEX
      (void) printf("  [%f+%fi", REAL(A->a[0]), IMAG(A->a[0]));
      for (j = 1; j < cols; j++)
	(void) printf(" %f+%fi", REAL(A->a[j*ldA]), IMAG(A->a[j*ldA]));
      (void) printf(" ;\n");
      for (i = 1; i < rows - 1; i++) {
	(void) printf("  %f+%fi", REAL(A->a[i]), IMAG(A->a[i]));
	for (j = 1; j < cols; j++)
	  (void) printf(" %f+%fi", REAL(A->a[i+j*ldA]), IMAG(A->a[i+j*ldA]));
	(void) printf(" ;\n");
      }
      (void) printf("  %f+%fi", REAL(A->a[i]), IMAG(A->a[i]));
      for (j = 1; j < cols; j++)
	(void) printf(" %f+%fi", REAL(A->a[i+j*ldA]), IMAG(A->a[i+j*ldA]));
      (void) printf("]\n");
#else
      (void) printf("  [%f", A->a[0]);
      for (j = 1; j < cols; j++)
	(void) printf(" %f", A->a[j*ldA]);
      (void) printf(" ;\n");
      for (i = 1; i < rows - 1; i++) {
	(void) printf("  %f", A->a[i]);
	for (j = 1; j < cols; j++)
	  (void) printf(" %f", A->a[i+j*ldA]);
	(void) printf(" ;\n");
      }
      (void) printf("  %f", A->a[i]);
      for (j = 1; j < cols; j++)
	(void) printf(" %f", A->a[i+j*ldA]);
      (void) printf("]\n");
#endif
    }
  }
}

real
check_ortho_amatrix(bool transA, pcamatrix A)
{
  pamatrix  B;
  uint      rows = A->rows;
  uint      cols = A->cols;
  real      norm, error;

  error = 0.0;

  if(transA) {
    B = new_amatrix(rows, rows);
    identity_amatrix(B);
    addmul_amatrix(-1.0, false, A, true, A, B);
    norm = normfrob_amatrix(B);
    del_amatrix(B);

    error = REAL_MAX(error, norm);
  }
  else {
    B = new_amatrix(cols, cols);
    identity_amatrix(B);
    addmul_amatrix(-1.0, true, A, false, A, B);
    norm = normfrob_amatrix(B);
    del_amatrix(B);

    error = REAL_MAX(error, norm);
  }

  return error;
}

/* ------------------------------------------------------------
 * Basic linear algebra
 * ------------------------------------------------------------ */

void
scale_amatrix(field alpha, pamatrix A)
{
  uint      rows = A->rows;
  uint      cols = A->cols;
  longindex ldA = A->ld;
  uint      j;

  for (j = 0; j < cols; j++)
    scal(rows, alpha, A->a + j *ldA, 1);
}

field
dotprod_amatrix(pcamatrix A, pcamatrix B)
{
  field     sum;
  uint      rows = A->rows;
  uint      cols = A->cols;
  longindex ldA = A->ld;
  longindex ldB = B->ld;
  uint      j;

  assert(rows == B->rows);
  assert(cols == B->cols);

  sum = 0.0;
  for (j = 0; j < cols; j++)
    sum += dot(true, false, rows, A->a + j * ldA, 1, B->a + j * ldB, 1);

  return sum;
}

real
normfrob_amatrix(pcamatrix A)
{
  real      sum;
  longindex ldA = A->ld;
  uint      j;

  sum = 0.0;
  for (j = 0; j < A->cols; j++)
    sum += REAL_SQR(nrm2(A->rows, A->a + j * ldA, 1));

  return REAL_SQRT(sum);
}

void
addeval_amatrix(field alpha, pcamatrix A, pcavector x, pavector y)
{
  assert(x->size == A->cols);
  assert(y->size == A->rows);

  if(A->rows == 0 || A->cols == 0)
    return;

  gemv(false, A->rows, A->cols, alpha, A->a, A->ld,
       x->v, 1, 1.0, y->v, 1);
}

void
addevaltrans_amatrix(field alpha, pcamatrix A, pcavector x, pavector y)
{
  assert(x->size == A->rows);
  assert(y->size == A->cols);

  if(A->rows == 0 || A->cols == 0)
    return;

  gemv(true, A->rows, A->cols, alpha, A->a, A->ld,
       x->v, 1, 1.0, y->v, 1);
}

void
mvm_amatrix(field alpha, bool transA, pcamatrix A,
	    pcavector x, pavector y)
{
  if(transA)
    addevaltrans_amatrix(alpha, A, x, y);
  else
    addeval_amatrix(alpha, A, x, y);
}

void
add_amatrix(field alpha, bool transA, pcamatrix A, pamatrix B)
{
  uint rows = B->rows;
  uint cols = B->cols;
  longindex ldA = A->ld;
  longindex ldB = B->ld;
  const field one = 1.0;
  int j;

  if(transA) {
    assert(rows == A->cols);
    assert(cols == A->rows);

    for(j=0; j<cols; j++)
      gemv(true, 1, rows, alpha, A->a+j, ldA, &one, 1, 1.0, B->a+j*ldB, 1);
  }
  else {
    assert(rows == A->rows);
    assert(cols == A->cols);

    for(j=0; j<cols; j++)
      axpy(rows, alpha, A->a+j*ldA, 1, B->a+j*ldB, 1);
  }
}

void
addmul_amatrix(field alpha, bool transA, pcamatrix A,
	       bool transB, pcamatrix B, pamatrix C)
{
  if (A->rows == 0 || A->cols == 0 || B->rows == 0 || B->cols == 0)
    return;

  if (transA) {
    if (transB) {
      /* C += alpha A^* B^* */
      assert(A->cols == C->rows);
      assert(B->rows == C->cols);
      assert(A->rows == B->cols);

      gemm(transA, transB, A->cols, B->rows, A->rows,
	   alpha, A->a, A->ld, B->a, B->ld, 1.0, C->a, C->ld);
    }
    else {
      /* C += alpha A^* B */
      assert(A->cols == C->rows);
      assert(B->cols == C->cols);
      assert(A->rows == B->rows);

      gemm(transA, transB, A->cols, B->cols, A->rows,
	   alpha, A->a, A->ld, B->a, B->ld, 1.0, C->a, C->ld);
    }
  }
  else {
    if (transB) {
      /* C += alpha A B^* */
      assert(A->rows == C->rows);
      assert(B->rows == C->cols);
      assert(A->cols == B->cols);

      gemm(transA, transB, A->rows, B->rows, A->cols,
	   alpha, A->a, A->ld, B->a, B->ld, 1.0, C->a, C->ld);
    }
    else {
      /* C += alpha A B */
      assert(A->rows == C->rows);
      assert(B->cols == C->cols);
      assert(A->cols == B->rows);

      gemm(transA, transB, A->rows, B->cols, A->cols,
	   alpha, A->a, A->ld, B->a, B->ld, 1.0, C->a, C->ld);
    }
  }
}

real
norm2_matrix(mvm_t mvm, void *A, uint rows, uint cols)
{
  avector   tmp1, tmp2;
  pavector  x, y;
  real      norm;
  uint      i;

  x = init_avector(&tmp1, cols);
  y = init_avector(&tmp2, rows);

  random_avector(x);
  norm = norm2_avector(x);
  i = 0;
  while (i < NORM_STEPS && norm > 0.0) {
    scale_avector(1.0 / norm, x);

    clear_avector(y);
    mvm(1.0, false, A, x, y);

    clear_avector(x);
    mvm(1.0, true, A, y, x);

    norm = norm2_avector(x);
    i++;
  }

  uninit_avector(y);
  uninit_avector(x);

  return REAL_SQRT(norm);
}

real
norm2_amatrix(pcamatrix A)
{
  return norm2_matrix((mvm_t) mvm_amatrix, (void *) A, A->rows,
		      A->cols);
}

real
norm2diff_matrix(mvm_t mvmA, void *A, mvm_t mvmB, void *B, uint rows,
		 uint cols)
{
  avector   tmp1, tmp2;
  pavector  x, y;
  real      norm;
  uint      i;

  x = init_avector(&tmp1, cols);
  y = init_avector(&tmp2, rows);

  random_avector(x);
  norm = norm2_avector(x);
  i = 0;
  while (i < NORM_STEPS && norm > 0.0) {
    scale_avector(1.0 / norm, x);

    clear_avector(y);
    mvmA(1.0, false, A, x, y);
    mvmB(-1.0, false, B, x, y);

    clear_avector(x);
    mvmA(1.0, true, A, y, x);
    mvmB(-1.0, true, B, y, x);

    norm = norm2_avector(x);
    i++;
  }

  uninit_avector(y);
  uninit_avector(x);

  return REAL_SQRT(norm);
}

real
norm2diff_amatrix(pcamatrix a, pcamatrix b)
{
  return norm2diff_matrix((mvm_t) mvm_amatrix, (void *) a,
			  (mvm_t) mvm_amatrix, (void *) b, a->rows,
			  a->cols);
}
