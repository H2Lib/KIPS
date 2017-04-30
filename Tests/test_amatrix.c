
#include <stdio.h>

#include "amatrix.h"
#include "settings.h"

#ifdef USE_FLOAT
static const real tolerance = 5.0e-5;
#else
static const real tolerance = 1.0e-12;
#endif

int
main()
{
  pamatrix  A, B, C, D;
  pavector  x, y, z;
  field     alpha;
  real      norm, error;
  uint      problems = 0;
  uint      rows, cols;
  uint      i, j, k;

  rows = 8;
  cols = 7;

  (void) printf("----------------------------------------\n"
		"Testing scale_amatrix\n");
  A = new_amatrix(rows, cols);
  random_amatrix(A);
  B = clone_amatrix(A);
  alpha = FIELD_RAND();
  scale_amatrix(alpha, B);
  error = norm = 0.0;
  for(j=0; j<cols; j++)
    for(i=0; i<rows; i++) {
      norm = REAL_MAX(norm, ABS(alpha * A->a[i+j*A->ld]));
      error = REAL_MAX(error, ABS(alpha * A->a[i+j*A->ld] - B->a[i+j*B->ld]));
    }
  (void) printf("  Error %.3e (%.3e)",
		error, (norm == 0.0 ? 0.0 : error/norm));
  if(error <= tolerance * norm)
    (void) printf("  --  Okay\n");
  else {
    (void) printf("  --  NOT Okay\a\n");
    problems++;
  }
  del_amatrix(B);
  del_amatrix(A);

  (void) printf("----------------------------------------\n"
		"Testing dotprod_amatrix\n");
  A = new_amatrix(rows, cols);
  random_amatrix(A);
  B = new_amatrix(rows, cols);
  random_amatrix(B);
  alpha = 0.0;
  for(j=0; j<cols; j++)
    for(i=0; i<rows; i++)
      alpha += CONJ(A->a[i+j*A->ld]) * B->a[i+j*B->ld];
  norm = ABS(alpha);
  alpha -= dotprod_amatrix(A, B);
  error = ABS(alpha);
  (void) printf("  Error %.3e (%.3e)",
		error, (norm == 0.0 ? 0.0 : error / norm));
  if(error <= tolerance * norm)
    (void) printf("  --  Okay\n");
  else {
    (void) printf("  --  NOT Okay\a\n");
    problems++;
  }
  del_amatrix(B);
  del_amatrix(A);

  (void) printf("----------------------------------------\n"
		"Testing addeval_amatrix\n");
  A = new_amatrix(rows, cols);
  random_amatrix(A);
  x = new_avector(cols);
  random_avector(x);
  y = new_avector(rows);
  random_avector(y);
  z = clone_avector(y);
  alpha = FIELD_RAND();
  for(j=0; j<cols; j++)
    for(i=0; i<rows; i++)
      y->v[i] += alpha * A->a[i+j*A->ld] * x->v[j];
  norm = norm2_avector(y);
  addeval_amatrix(-alpha, A, x, y);
  add_avector(-1.0, y, z);
  error = norm2_avector(z);
  (void) printf("  Error %.3e (%.3e)",
		error, (norm == 0.0 ? 0.0 : error / norm));
  if(error <= tolerance * norm)
    (void) printf("  --  Okay\n");
  else {
    (void) printf("  --  NOT Okay\a\n");
    problems++;
  }
  del_avector(z);
  del_avector(y);
  del_avector(x);
  del_amatrix(A);

  (void) printf("----------------------------------------\n"
		"Testing addevaltrans_amatrix\n");
  A = new_amatrix(rows, cols);
  random_amatrix(A);
  x = new_avector(rows);
  random_avector(x);
  y = new_avector(cols);
  random_avector(y);
  z = clone_avector(y);
  alpha = FIELD_RAND();
  for(i=0; i<cols; i++)
    for(j=0; j<rows; j++)
      y->v[i] += alpha * CONJ(A->a[j+i*A->ld]) * x->v[j];
  norm = norm2_avector(y);
  addevaltrans_amatrix(-alpha, A, x, y);
  add_avector(-1.0, y, z);
  error = norm2_avector(z);
  (void) printf("  Error %.3e (%.3e)",
		error, (norm == 0.0 ? 0.0 : error / norm));
  if(error <= tolerance * norm)
    (void) printf("  --  Okay\n");
  else {
    (void) printf("  --  NOT Okay\a\n");
    problems++;
  }
  del_avector(z);
  del_avector(y);
  del_avector(x);
  del_amatrix(A);

  (void) printf("----------------------------------------\n"
		"Testing add_amatrix\n");
  A = new_amatrix(rows, cols);
  random_amatrix(A);
  B = new_amatrix(rows, cols);
  random_amatrix(B);
  C = clone_amatrix(A);
  alpha = FIELD_RAND();
  for(j=0; j<cols; j++)
    for(i=0; i<rows; i++)
      C->a[i+j*C->ld] += alpha * B->a[i+j*B->ld];
  norm = normfrob_amatrix(C);
  add_amatrix(-alpha, false, B, C);
  add_amatrix(-1.0, false, A, C);
  error = normfrob_amatrix(C);
  (void) printf("  Error %.3e (%.3e)",
		error, (norm == 0.0 ? 0.0 : error / norm));
  if(error <= tolerance * norm)
    (void) printf("  --  Okay\n");
  else {
    (void) printf("  --  NOT Okay\a\n");
    problems++;
  }
  del_amatrix(C);
  del_amatrix(B);
  del_amatrix(A);

  A = new_amatrix(rows, cols);
  random_amatrix(A);
  B = new_amatrix(cols, rows);
  random_amatrix(B);
  C = clone_amatrix(A);
  alpha = FIELD_RAND();
  for(j=0; j<cols; j++)
    for(i=0; i<rows; i++)
      C->a[i+j*C->ld] += alpha * CONJ(B->a[j+i*B->ld]);
  norm = normfrob_amatrix(C);
  add_amatrix(-alpha, true, B, C);
  add_amatrix(-1.0, false, A, C);
  error = normfrob_amatrix(C);
  (void) printf("  Error %.3e (%.3e)",
		error, (norm == 0.0 ? 0.0 : error / norm));
  if(error <= tolerance * norm)
    (void) printf("  --  Okay\n");
  else {
    (void) printf("  --  NOT Okay\a\n");
    problems++;
  }
  del_amatrix(C);
  del_amatrix(B);
  del_amatrix(A);

  (void) printf("----------------------------------------\n"
		"Testing addmul_amatrix\n");
  A = new_amatrix(rows, cols);
  random_amatrix(A);
  B = new_amatrix(cols, rows);
  random_amatrix(B);
  C = new_amatrix(rows, rows);
  random_amatrix(C);
  D = clone_amatrix(C);
  alpha = FIELD_RAND();
  for(j=0; j<D->cols; j++)
    for(i=0; i<D->rows; i++)
      for(k=0; k<A->cols; k++)
	D->a[i+j*D->ld] += alpha * A->a[i+k*A->ld] * B->a[k+j*B->ld];
  norm = normfrob_amatrix(D);
  addmul_amatrix(-alpha, false, A, false, B, D);
  add_amatrix(-1.0, false, C, D);
  error = normfrob_amatrix(D);
  (void) printf("  Error %.3e (%.3e)",
		error, (norm == 0.0 ? 0.0 : error / norm));
  if(error <= tolerance * norm)
    (void) printf("  --  Okay\n");
  else {
    (void) printf("  --  NOT Okay\a\n");
    problems++;
  }
  del_amatrix(D);
  del_amatrix(C);
  del_amatrix(B);
  del_amatrix(A);

  A = new_amatrix(cols, rows);
  random_amatrix(A);
  B = new_amatrix(cols, rows);
  random_amatrix(B);
  C = new_amatrix(rows, rows);
  random_amatrix(C);
  D = clone_amatrix(C);
  alpha = FIELD_RAND();
  for(j=0; j<D->cols; j++)
    for(i=0; i<D->rows; i++)
      for(k=0; k<A->rows; k++)
	D->a[i+j*D->ld] += alpha * CONJ(A->a[k+i*A->ld]) * B->a[k+j*B->ld];
  norm = normfrob_amatrix(D);
  addmul_amatrix(-alpha, true, A, false, B, D);
  add_amatrix(-1.0, false, C, D);
  error = normfrob_amatrix(D);
  (void) printf("  Error %.3e (%.3e)",
		error, (norm == 0.0 ? 0.0 : error / norm));
  if(error <= tolerance * norm)
    (void) printf("  --  Okay\n");
  else {
    (void) printf("  --  NOT Okay\a\n");
    problems++;
  }
  del_amatrix(D);
  del_amatrix(C);
  del_amatrix(B);
  del_amatrix(A);

  A = new_amatrix(rows, cols);
  random_amatrix(A);
  B = new_amatrix(rows, cols);
  random_amatrix(B);
  C = new_amatrix(rows, rows);
  random_amatrix(C);
  D = clone_amatrix(C);
  alpha = FIELD_RAND();
  for(j=0; j<D->cols; j++)
    for(i=0; i<D->rows; i++)
      for(k=0; k<A->cols; k++)
	D->a[i+j*D->ld] += alpha * A->a[i+k*A->ld] * CONJ(B->a[j+k*B->ld]);
  norm = normfrob_amatrix(D);
  addmul_amatrix(-alpha, false, A, true, B, D);
  add_amatrix(-1.0, false, C, D);
  error = normfrob_amatrix(D);
  (void) printf("  Error %.3e (%.3e)",
		error, (norm == 0.0 ? 0.0 : error / norm));
  if(error <= tolerance * norm)
    (void) printf("  --  Okay\n");
  else {
    (void) printf("  --  NOT Okay\a\n");
    problems++;
  }
  del_amatrix(D);
  del_amatrix(C);
  del_amatrix(B);
  del_amatrix(A);

  A = new_amatrix(cols, rows);
  random_amatrix(A);
  B = new_amatrix(rows, cols);
  random_amatrix(B);
  C = new_amatrix(rows, rows);
  random_amatrix(C);
  D = clone_amatrix(C);
  alpha = FIELD_RAND();
  for(j=0; j<D->cols; j++)
    for(i=0; i<D->rows; i++)
      for(k=0; k<A->rows; k++)
	D->a[i+j*D->ld] += alpha * CONJ(A->a[k+i*A->ld]) * CONJ(B->a[j+k*B->ld]);
  norm = normfrob_amatrix(D);
  addmul_amatrix(-alpha, true, A, true, B, D);
  add_amatrix(-1.0, false, C, D);
  error = normfrob_amatrix(D);
  (void) printf("  Error %.3e (%.3e)",
		error, (norm == 0.0 ? 0.0 : error / norm));
  if(error <= tolerance * norm)
    (void) printf("  --  Okay\n");
  else {
    (void) printf("  --  NOT Okay\a\n");
    problems++;
  }
  del_amatrix(D);
  del_amatrix(C);
  del_amatrix(B);
  del_amatrix(A);
  (void) printf("----------------------------------------\n"
		"  %u matrices and\n"
		"  %u vectors still active\n"
		"  %u errors found\n", getactives_amatrix(),
		getactives_avector(), problems);

  return problems;
}
