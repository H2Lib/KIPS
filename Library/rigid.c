#include "rigid.h"
#include "blas.h"
#include <stdlib.h>

/** ---------------------------------------------------------------------
 *    Constructors and destructors
 *  --------------------------------------------------------------------- */

pmolecule
new_molecule (uint n) {
  pmolecule mol = (pmolecule) allocmem (sizeof (molecule));
  mol->n = n;
  mol->R = new_amatrix (3, 3);
  mol->q = new_quaternion ();
  mol->qv = new_quaternion ();
  mol->qa = new_quaternion ();
  mol->r = (pavector *) allocmem (sizeof (pavector) * n);
  mol->rp = (pavector *) allocmem (sizeof (pavector) * n);
  for (uint i=0; i<n; i++) {
    mol->r[i] = new_avector (3);
    mol->rp[i] = new_avector (3);
  }
  mol->v = new_avector (3);
  mol->f = new_avector (3);
  mol->w = allocreal (3);
  mol->t = new_avector (3);
  mol->tp = new_avector (3);
  mol->com = new_avector (3);
  mol->I = allocreal (3);
  
  
  return mol;
}

void
del_molecule (pmolecule mol) {
  del_amatrix (mol->R);
  del_quaternion (mol->q);
  del_quaternion (mol->qv);
  del_quaternion (mol->qa);
  freemem (mol->I);
  freemem (mol->w);
  del_avector (mol->t);
  del_avector (mol->tp);
  del_avector (mol->v);
  del_avector (mol->f);
  del_avector (mol->com);
  for (uint i=0; i<mol->n; i++) {
    del_avector (mol->r[i]);
    del_avector (mol->rp[i]);
  }
  freemem (mol->r);
  freemem (mol->rp);
  freemem (mol);
}


/** ---------------------------------------------------------------------
 *    Reference orientation
 *  --------------------------------------------------------------------- */

void
refOrientation_molecule (pmolecule mol) {
  pamatrix R = mol->R;
  pquaternion q = mol->q;
  preal I = mol->I;
  pfield a = R->a;
  
  // Diagonalization
  Jacobi_quaternion (R, q);
  
  // Principal moments of inertia
  I[0] = a[0];
  I[1] = a[4];
  I[2] = a[8];
  assert (REAL_ABS (I[0]) > KIPS_ALMOST_ZERO);
  assert (REAL_ABS (I[1]) > KIPS_ALMOST_ZERO);
  assert (REAL_ABS (I[2]) > KIPS_ALMOST_ZERO);
  
  // Rotation matrix
  buildRotation_quaternion (q, R);
  
  // Reference frame
  clear_avector (mol->tp);
  addevaltrans_amatrix (1.0, R, mol->t, mol->tp);
  for (uint i=0; i<mol->n; i++) {
    clear_avector (mol->rp[i]);
    addevaltrans_amatrix (1.0, R, mol->r[i], mol->rp[i]);
  }
}

void
adjust_molecule (pcquaternion q, pmolecule mol) {
  uint i, n;
  pamatrix R = mol->R;
  pavector *r, *rp;
  r = mol->r;
  rp = mol->rp;
  n = mol->n;
  
  buildRotation_quaternion (q, R);
  for (i=0; i<n; i++) {
    clear_avector (r[i]);
    addeval_amatrix (1.0, R, rp[i], r[i]);
  }
}

/** ---------------------------------------------------------------------
 *    Time-stepping methods
 *  --------------------------------------------------------------------- */

void
angularVelocity_molecule (pmolecule mol) {
  pquaternion q, qv;
  preal w = mol->w;
  q = mol->q;
  qv = mol->qv;
  
  w[0] = 2.0 * (-q[1] * qv[0] + q[0] * qv[1] + q[3] * qv[2] - q[2] * qv[3]);
  w[1] = 2.0 * (-q[2] * qv[0] - q[3] * qv[1] + q[0] * qv[2] + q[1] * qv[3]);
  w[2] = 2.0 * (-q[3] * qv[0] + q[2] * qv[1] - q[1] * qv[2] + q[0] * qv[3]);
}

void
quaternionAcceleration_molecule (pmolecule mol) {
  pquaternion q, qv, qa;
  preal w, I;
  real qv2, wvx, wvy, wvz;
  pfield tp = mol->tp->v;
  w = mol->w;
  I = mol->I;
  q = mol->q;
  qv = mol->qv;
  qa = mol->qa;
  
  qv2 = qv[0] * qv[0] + qv[1] * qv[1] + qv[2] * qv[2] + qv[3] * qv[3];
  wvx = tp[0] + w[1] * w[2] * (I[1] - I[2]) / I[0];
  wvy = tp[1] + w[2] * w[0] * (I[2] - I[0]) / I[2];
  wvz = tp[2] + w[0] * w[1] * (I[0] - I[1]) / I[2];
  qa[0] = -q[0] * qv2 - 0.5 * (q[1] * wvx + q[2] * wvy + q[3] * wvz);
  qa[1] = -q[1] * qv2 + 0.5 * (q[0] * wvx - q[3] * wvy + q[2] * wvz);
  qa[2] = -q[2] * qv2 + 0.5 * (q[3] * wvx + q[0] * wvy - q[1] * wvz);
  qa[3] = -q[3] * qv2 + 0.5 * (-q[2] * wvx + q[1] * wvy + q[0] * wvz);
}
