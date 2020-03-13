#include "quaternion.h"
#include "basic.h"
#include <assert.h>


/* ---------------------------------------------------------------------
 *    Constructors and destructors
 * --------------------------------------------------------------------- */

pquaternion 
new_quaternion () {
  return allocreal (4);
}

void
del_quaternion (pquaternion q) {
  freemem (q);
}


/* ---------------------------------------------------------------------
 *    Basic operations
 * --------------------------------------------------------------------- */

void 
clear_quaternion (pquaternion q) {
  q[0] = 0.0;
  q[1] = 0.0;
  q[2] = 0.0;
  q[3] = 0.0;
}

void
copy_quaternion (pcquaternion p, pquaternion q) {
  q[0] = p[0];
  q[1] = p[1];
  q[2] = p[2];
  q[3] = p[3];
}

void
normalize_quaternion (pquaternion q) {
  real norm;
  
  norm = nrm2 (4, q, 1);
  assert (norm > KIPS_ALMOST_ZERO);
  norm = 1.0 / norm;
  scal (4, norm, q, 1);
}

void update_quaternion (pcquaternion p, pquaternion q) {
  real z0, z1, z2, z3;
  z0 = q[0] * p[0] - q[1] * p[1] - q[2] * p[2] - q[3] * p[3];
  z1 = q[0] * p[1] + q[1] * p[0] + q[2] * p[3] - q[3] * p[2];
  z2 = q[0] * p[2] - q[1] * p[3] + q[2] * p[0] + q[3] * p[1];
  z3 = q[0] * p[3] + q[1] * p[2] - q[2] * p[1] + q[3] * p[0];
  
  q[0] = z0;
  q[1] = z1;
  q[2] = z2;
  q[3] = z3;
}

void buildRotation_quaternion (pcquaternion q, preal R) {
  
  R[0] = q[0] * q[0] + q[1] * q[1] - q[2] * q[2] - q[3] * q[3];
  R[1] = 2.0 * (q[1] * q[2] + q[0] * q[3]);
  R[2] = 2.0 * (q[1] * q[3] - q[0] * q[2]);
  R[3] = 2.0 * (q[1] * q[2] - q[0] * q[3]);
  R[4] = q[0] * q[0] - q[1] * q[1] + q[2] * q[2] - q[3] * q[3];
  R[5] = 2.0 * (q[2] * q[3] + q[0] * q[1]);
  R[6] = 2.0 * (q[1] * q[3] + q[0] * q[2]);
  R[7] = 2.0 * (q[2] * q[3] - q[0] * q[1]);
  R[8] = q[0] * q[0] - q[1] * q[1] - q[2] * q[2] + q[3] * q[3];
}


/* ---------------------------------------------------------------------
 *    Diagonalization
 * --------------------------------------------------------------------- */

uint Jacobi_quaternion (preal A, pquaternion q) {
  
  uint i, j, k, l, iter;
  real tau, t, c, s, h, off, b, c_h, s_h, t_h, normf;
  preal p;
  
  p = new_quaternion ();
  q[0] = 1.0;
  q[1] = 0.0;
  q[2] = 0.0;
  q[3] = 0.0;
  iter = 0;
  
  // Determine squared Frobenius norm of the off-diagonal part.
  off = 0.0;
  normf = 0.0;
  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      b = ABSSQR (A[i+j*3]);
      normf += b;
      if (i != j) {
        off += b;
      }
    }
  }
  
  while (off/normf > KIPS_MACH_EPS) {
    // Determine largest off-diagonal matrix entry.
    i = 0;
    j = 1;
    for (l=1; l<3; l++) {
      for (k=0; k<l; k++) {
        if (ABS (A[k+l*3]) > ABS (A[i+j*3])) {
          i = k;
          j = l;
        }
      }
    }
    
    // Determine suitable Givens rotation.
    b = A[i+j*3];
    tau = (A[j+j*3] - A[i+i*3]) / (2.0 * b);
    t = -SIGN1 (tau) / (ABS (tau) + REAL_SQRT (ABSSQR (tau) + 1));
    c = 1.0 / REAL_SQRT (1 + ABSSQR (t));
    s = t * c;
    
    // Half angle formulas.
    t_h = t / (1.0 + REAL_SQRT (1.0 + ABSSQR (t)));
    c_h = REAL_SQRT ((1.0 + c) / 2.0);
    s_h = t_h * c_h;
    
    // Quaternion corresponding to the current Givens rotation.
    p[0] = c_h;
    p[i+1] = 0.0;
    p[j+1] = 0.0;
    k = 0;
    if (k == i) {
      k = 1;
      if (k == j) {
        k = 2;
      }
    }
    p[k+1] = (k == 1) ? -s_h : s_h;
    
    // Multiply quaternions and update matrix.
    update_quaternion (p, q);
    for (k=0; k<3; k++) {
      h = A[k+i*3];
      A[k+i*3] = h * c + A[k+j*3] * s;
      A[k+j*3] = A[k+j*3] * c - h * CONJ (s);
    }
    for (k=0; k<3; k++) {
      h = A[i+k*3];
      A[i+k*3] = h * c + A[j+k*3] * CONJ (s);
      A[j+k*3] = A[j+k*3] * c - h * s;
    }
    off -= 2.0 * ABSSQR (b);
    iter++;
  }
  
  del_quaternion (p);
  
  return iter;
}
