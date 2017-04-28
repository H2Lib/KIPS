
#include <stdio.h>

#include "amatrix.h"
#include "settings.h"

#ifdef USE_FLOAT
static const real tolerance = 5.0e-5;
#else
static const real tolerance = 1.0e-12;
#endif

#define IS_IN_RANGE(a, b, c) (((a) <= (b)) && ((b) <= (c)))

#ifdef USE_COMPLEX
static field alpha = 1.0 + 1.0 * I;
#else
static field alpha = 1.0;
#endif

int
main()
{
  pamatrix  a, acopy, l, ld, dlt, ldltcopy, r, q, qr, X, B;
  pavector  x, b, tau, lvec, dvec;
  uint     problems = 0;
  field    *data;
  uint      rows, cols, rhs, i, j;
  uint     *colpiv;
  real      error;

  rows = 8;
  cols = 7;

  (void) printf("----------------------------------------\n"
		"Call utility functions\n");
  a = new_amatrix(rows, cols);
  random_amatrix(a);

  (void) printf("Print matrix to stdout\n");
  print_amatrix(a);
  (void) printf("Print matrix in matlab format\n");
  print_matlab_amatrix(a);

  del_amatrix(a);

  a = new_amatrix(rows, 0);
  random_amatrix(a);

  (void) printf("Print matrix to stdout\n");
  print_amatrix(a);
  (void) printf("Print matrix in matlab format\n");
  print_matlab_amatrix(a);

  del_amatrix(a);

  a = new_amatrix(0, cols);
  random_amatrix(a);

  (void) printf("Print matrix to stdout\n");
  print_amatrix(a);
  (void) printf("Print matrix in matlab format\n");
  print_matlab_amatrix(a);

  del_amatrix(a);

  a = new_amatrix(1, 0);
  random_amatrix(a);

  (void) printf("Print matrix to stdout\n");
  print_amatrix(a);
  (void) printf("Print matrix in matlab format\n");
  print_matlab_amatrix(a);

  del_amatrix(a);

  a = new_amatrix(1, cols);
  random_amatrix(a);

  (void) printf("Print matrix to stdout\n");
  print_amatrix(a);
  (void) printf("Print matrix in matlab format\n");
  print_matlab_amatrix(a);

  del_amatrix(a);

  (void) printf("----------------------------------------\n"
		"  %u matrices and\n"
		"  %u vectors still active\n"
		"  %u errors found\n", getactives_amatrix(),
		getactives_avector(), problems);

  return problems;
}
