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

void buildRotation_quaternion (pcquaternion q, pamatrix R) {
  assert (R->rows == 3 && R->cols == 3);
  preal r = R->a;
  
  r[0] = q[0] * q[0] + q[1] * q[1] - q[2] * q[2] - q[3] * q[3];
  r[1] = 2.0 * (q[1] * q[2] + q[0] * q[3]);
  r[2] = 2.0 * (q[1] * q[3] - q[0] * q[2]);
  r[3] = 2.0 * (q[1] * q[2] - q[0] * q[3]);
  r[4] = q[0] * q[0] - q[1] * q[1] + q[2] * q[2] - q[3] * q[3];
  r[5] = 2.0 * (q[2] * q[3] + q[0] * q[1]);
  r[6] = 2.0 * (q[1] * q[3] + q[0] * q[2]);
  r[7] = 2.0 * (q[2] * q[3] - q[0] * q[1]);
  r[8] = q[0] * q[0] - q[1] * q[1] - q[2] * q[2] + q[3] * q[3];
}


/* ---------------------------------------------------------------------
 *    Diagonalization
 * --------------------------------------------------------------------- */

uint Jacobi_quaternion (pamatrix A, pquaternion q) {
  
  uint lda, n, i, j, k, l, iter;
  real tau, t, c, s, h, off, b, c_h, s_h, t_h, normf;
  pfield aa;
  preal p;
  
  n = 3;
  lda = A->ld;
  aa = A->a;
  
  assert (n == A->cols);
  assert (n == A->rows);
  
  p = new_quaternion ();
  q[0] = 1.0;
  q[1] = 0.0;
  q[2] = 0.0;
  q[3] = 0.0;
  iter = 0;
  
  // Determine squared Frobenius norm of the off-diagonal part.
  off = 0.0;
  normf = 0.0;
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) {
      b = ABSSQR (aa[i+j*lda]);
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
    for (l=1; l<n; l++) {
      for (k=0; k<l; k++) {
        if (ABS (aa[k+l*lda]) > ABS (aa[i+j*lda])) {
          i = k;
          j = l;
        }
      }
    }
    
    // Determine suitable Givens rotation.
    b = aa[i+j*lda];
    tau = (aa[j+j*lda] - aa[i+i*lda]) / (2.0 * b);
    t = -SIGN1 (tau) / (ABS (tau) + REAL_SQRT (ABSSQR (tau) + 1));
    c = 1.0 / REAL_SQRT (1 + t * CONJ (t));
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
    for (k=0; k<n; k++) {
      h = aa[k+i*lda];
      aa[k+i*lda] = h * c + aa[k+j*lda] * s;
      aa[k+j*lda] = aa[k+j*lda] * c - h * CONJ (s);
    }
    for (k=0; k<n; k++) {
      h = aa[i+k*lda];
      aa[i+k*lda] = h * c + aa[j+k*lda] * CONJ (s);
      aa[j+k*lda] = aa[j+k*lda] * c - h * s;
    }
    off -= 2.0 * ABSSQR (b);
    iter++;
  }
  
  del_quaternion (p);
  
  return iter;
}
