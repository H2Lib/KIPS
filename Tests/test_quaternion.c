#include "amatrix.h"
#include "quaternion.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>


int main (void) {
  pamatrix A, R;
  preal a;
  pquaternion q;
  uint iter;
  real e1, e2, e3;
  
  A = new_amatrix (3, 3);
  R = new_amatrix (3, 3);
  q = new_quaternion ();
  a = A->a;
  
  printf ("Testing quaternion rotation by diagonalization of a real symmetric 3x3-matrix.\n");
  
  a[0] = 1.0; a[1] = 2.0; a[2] = 0.0;
  a[3] = a[1]; a[4] = 1.0; a[5] = 2.0;
  a[6] = a[2]; a[7] = a[5]; a[8] = 1.0;
  printf ("Non-diagonal matrix:\n");
  print_amatrix (A);
  
  iter = Jacobi_quaternion (A, q);
  printf ("Diagonalization complete after %u steps\n", iter);
  printf ("Diagonal matrix:\n");
  print_amatrix (A);
  
  buildRotation_quaternion (q, R);
  printf ("Corresponding rotation matrix:\n");
  print_amatrix (R);
  
  e1 = 1.0 - REAL_SQRT (8.0);
  e2 = 1.0 + REAL_SQRT (8.0);
  e3 = 1.0;
  
  e1 = REAL_ABS (a[0] - e1) / REAL_ABS (e1);
  e2 = REAL_ABS (a[4] - e2) / REAL_ABS (e2);
  e3 = REAL_ABS (a[8] - e3) / REAL_ABS (e3);
  
  printf ("Relative errors of the eigenvalues:\n");
  printf ("e_1 = %e   e_2 = %e   e_3 = %e\n", e1, e2, e3);
  
  del_quaternion (q);
  del_amatrix (A);
  del_amatrix (R);
  return EXIT_SUCCESS;
}
