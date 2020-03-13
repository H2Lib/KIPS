#include "amatrix.h"
#include "quaternion.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>


int main (void) {
  pamatrix A, R;
  preal a, r;
  pquaternion q;
  uint iter;
  
  A = new_amatrix (3, 3);
  R = new_amatrix (3, 3);
  q = new_quaternion ();
  a = A->a;
  r = R->a;
  
  a[0] = 1.0; a[1] = 4.0; a[2] = 5.0;
  a[3] = a[1]; a[4] = 2.0; a[5] = 6.0;
  a[6] = a[2]; a[7] = a[5]; a[8] = 3.0;
  printf ("Non-diagonal matrix:\n");
  print_amatrix (A);
  
  iter = Jacobi_quaternion (a, q);
  printf ("Diagonalization complete after %u steps\n", iter);
  printf ("Diagonal matrix:\n");
  print_amatrix (A);
  
  buildRotation_quaternion (q, r);
  printf ("Corresponding rotation matrix:\n");
  print_amatrix (R);
  
  del_quaternion (q);
  del_amatrix (A);
  del_amatrix (R);
  return EXIT_SUCCESS;
}
