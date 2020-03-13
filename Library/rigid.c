#include "rigid.h"
#include "blas.h"
#include <stdlib.h>
#include <stdio.h>

/** ---------------------------------------------------------------------
 *    Constructors and destructors
 *  --------------------------------------------------------------------- */

static void 
set_zero (uint n, preal x) {
  for (uint i=0; i<n; i++) 
    x[i] = 0.0;
}

pmolecule
new_molecule (uint n) {
  pmolecule mol = (pmolecule) allocmem (sizeof (molecule));
  mol->n = n;
  mol->R = allocreal (9);
  mol->q = new_quaternion ();
  mol->qv = new_quaternion ();
  mol->qa = new_quaternion ();
  mol->qh = new_quaternion ();
  mol->r = (preal *) allocmem (sizeof (preal) * n);
  mol->rp = (preal *) allocmem (sizeof (preal) * n);
  mol->r[0] = allocreal (3 * n);
  mol->rp[0] = allocreal (3 * n);
  for (uint i=1; i<n; i++) {
    mol->r[i] = mol->r[0] + 3 * i;
    mol->rp[i] = mol->rp[0] + 3 * i;
  }
  mol->v = allocreal (3);
  mol->f = allocreal (3);
  mol->w = allocreal (3);
  mol->t = allocreal (3);
  mol->tp = allocreal (3);
  mol->com = allocreal (3);
  mol->I = allocreal (3);
  
  return mol;
}

void
del_molecule (pmolecule mol) {
  freemem (mol->R);
  del_quaternion (mol->q);
  del_quaternion (mol->qv);
  del_quaternion (mol->qa);
  del_quaternion (mol->qh);
  freemem (mol->I);
  freemem (mol->w);
  freemem (mol->t);
  freemem (mol->tp);
  freemem (mol->v);
  freemem (mol->f);
  freemem (mol->com);
  freemem (mol->r[0]);
  freemem (mol->rp[0]);
  freemem (mol->r);
  freemem (mol->rp);
  freemem (mol);
}


/** ---------------------------------------------------------------------
 *    Reference orientation
 *  --------------------------------------------------------------------- */

void
initRefOrientationFromFile_molecule (const char *file, uint n, pmolecule *mol) {
  real Ix, Iy, Iz, qw, qi, qj, qk;
  preal I;
  pquaternion q;
  uint i, j, k, N;
  FILE *stream;
  char buf[100];
  
  stream = fopen (file, "r");
  if (stream == NULL) {
    printf ("Error: File not found.\n");
    exit (0);
  }
  fscanf (stream, "%u", &N);
  assert (n == N);
  fgets (buf, 100, stream);
  for (i=0; i<n; i++) {
    #ifdef USE_FLOAT
      fscanf (stream, "%u %u %e %e %e %e %e %e %e\n", &j, &k, &Ix, &Iy, &Iz, &qw, &qi, &qj, &qk);
    #else
      fscanf (stream, "%u %u %le %le %le %le %le %le %le\n", &j, &k, &Ix, &Iy, &Iz, &qw, &qi, &qj, &qk);
    #endif
    k--;
    I = mol[k]->I;
    q = mol[k]->q;
    I[0] = Ix;
    I[1] = Iy;
    I[2] = Iz;
    q[0] = qw;
    q[1] = qi;
    q[2] = qj;
    q[3] = qk;
    buildRotation_quaternion (q, mol[k]->R);
  }
}

void
refOrientation_molecule (pmolecule mol) {
  preal t, tp, R, I, *rp, *r;
  pquaternion q = mol->q;
  uint i, n;
  t = mol->t;
  tp = mol->tp;
  R = mol->R;
  I = mol->I;
  n = mol->n;
  r = mol->r;
  rp = mol->rp;
  
  // Diagonalization
  Jacobi_quaternion (R, q);
  
  // Principal moments of inertia
  I[0] = R[0];
  I[1] = R[4];
  I[2] = R[8];
  assert (REAL_ABS (I[0]) > KIPS_ALMOST_ZERO);
  assert (REAL_ABS (I[1]) > KIPS_ALMOST_ZERO);
  assert (REAL_ABS (I[2]) > KIPS_ALMOST_ZERO);
  
  // Rotation matrix
  buildRotation_quaternion (q, R);
  
  // Reference frame
  set_zero (3 * n, rp[0]);
  for (i=0; i<n; i++) {
    gemv (true, 3, 3, 1.0, R, 3, r[i], 1, 1.0, rp[i], 1);
  }
  set_zero (3, tp);
  gemv (true, 3, 3, 1.0, R, 3, t, 1, 1.0, tp, 1);
}

void
adjust_molecule (pcquaternion q, pmolecule mol) {
  uint i, n;
  preal R, *r, *rp;
  r = mol->r;
  rp = mol->rp;
  n = mol->n;
  R = mol->R;
  
  buildRotation_quaternion (q, R);
  set_zero (3 * n, r[0]);
  for (i=0; i<n; i++) {
    gemv (false, 3, 3, 1.0, R, 3, rp[i], 1, 1.0, r[i], 1);
  }
}

/** ---------------------------------------------------------------------
 *    Time-stepping methods
 *  --------------------------------------------------------------------- */

static void
calc_quatvel_from_angvel (pcreal w, pcquaternion q, pquaternion qv) {
  qv[0] = -0.5 * (q[1] * w[0] + q[2] * w[1] + q[3] * w[2]);
  qv[1] = 0.5 * (q[0] * w[0] - q[3] * w[1] + q[2] * w[2]);
  qv[2] = 0.5 * (q[3] * w[0] + q[0] * w[1] - q[1] * w[2]);
  qv[3] = 0.5 * (-q[2] * w[0] + q[1] * w[1] + q[0] * w[2]); 
}

real 
initVelocitiesFromFile_molecule (const char *file, uint n, pmolecule *mol) {
  real vx, vy, vz, E, E_rot, E_trans, m;
  preal I, v, w, wp, R;
  uint i, j, k, N;
  FILE *stream;
  char buf[100];
  
  E = 0.0;
  E_rot = 0.0;
  E_trans = 0.0;
  w = allocreal (3);
  
  stream = fopen (file, "r");
  if (stream == NULL) {
    printf ("Error: File not found.\n");
    exit (0);
  }
  fscanf (stream, "%u", &N);
  assert (n == N);
  fgets (buf, 100, stream);
  for (i=0; i<n; i++) {
    #ifdef USE_FLOAT
      fscanf (stream, "%u %u %e %e %e %e %e %e\n", &j, &k, &vx, &vy, &vz, &w, &wy, &wz);
    #else
      fscanf (stream, "%u %u %le %le %le %le %le %le\n", &j, &k, &vx, &vy, &vz, &w[0], &w[1], &w[2]);
    #endif
    k--;
    v = mol[k]->v;
    wp = mol[k]->w;
    I = mol[k]->I;
    m = mol[k]->m;
    R = mol[k]->R;
    v[0] = vx; v[1] = vy; v[2] = vz;
    set_zero (3, wp);
    gemv (true, 3, 3, 1.0, R, 3, w, 1, 1.0, wp, 1);
    //printf ("%lf  %lf  %lf  %lf\n", mol[k]->q[0], mol[k]->q[1], mol[k]->q[2], mol[k]->q[3]);
    //printf ("%le %lf %lf %lf %le %le %le %le %le %le\n", m, v[0], v[1], v[2], I[0], I[1], I[2], w[0], w[1], w[2]);
    calc_quatvel_from_angvel (wp, mol[k]->q, mol[k]->qv);
    E_rot += 0.5 * (wp[0] * wp[0] * I[0] + wp[1] * wp[1] * I[1] 
            + wp[2] * wp[2] * I[2]);
    E_trans += 0.5 * m * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  }
  E = E_rot + E_trans;
  
  freemem (w);
  printf ("translational energy: %le  rotational energy: %le  total kinetic energy: %le\n", E_trans, E_rot, E);
  //printf ("Initial temperature: %f\n", 2.0 * E / ((6 * n - 3) * K_B));
  return E;
}

static void
angular_velocity_molecule (pcquaternion q, pcquaternion qv, preal w) {
  w[0] = 2.0 * (-q[1] * qv[0] + q[0] * qv[1] + q[3] * qv[2] - q[2] * qv[3]);
  w[1] = 2.0 * (-q[2] * qv[0] - q[3] * qv[1] + q[0] * qv[2] + q[1] * qv[3]);
  w[2] = 2.0 * (-q[3] * qv[0] + q[2] * qv[1] - q[1] * qv[2] + q[0] * qv[3]);
}

static void
quaternion_acceleration_molecule (pcreal w, pcreal I, pcreal tp, 
                                  pcquaternion q, pcquaternion qv, pquaternion qa) {
  real qv2, wvx, wvy, wvz;
  
  qv2 = qv[0] * qv[0] + qv[1] * qv[1] + qv[2] * qv[2] + qv[3] * qv[3];
  wvx = tp[0] + w[1] * w[2] * (I[1] - I[2]) / I[0];
  wvy = tp[1] + w[2] * w[0] * (I[2] - I[0]) / I[2];
  wvz = tp[2] + w[0] * w[1] * (I[0] - I[1]) / I[2];
  qa[0] = -q[0] * qv2 - 0.5 * (q[1] * wvx + q[2] * wvy + q[3] * wvz);
  qa[1] = -q[1] * qv2 + 0.5 * (q[0] * wvx - q[3] * wvy + q[2] * wvz);
  qa[2] = -q[2] * qv2 + 0.5 * (q[3] * wvx + q[0] * wvy - q[1] * wvz);
  qa[3] = -q[3] * qv2 + 0.5 * (-q[2] * wvx + q[1] * wvy + q[0] * wvz);
}

void 
positionVerlet_molecule (real delta, pmolecule mol) {
  pquaternion q, qv, qa, qh;
  real m = mol->m;
  preal I, f, v, w, com, tp;
  
  w = mol->w;
  I = mol->I;
  f = mol->f;
  v = mol->v;
  com = mol->com;
  tp = mol->tp;
  q = mol->q;
  qv = mol->qv;
  qa = mol->qa;
  qh = mol->qh;
  
  // Translational motion
  
  // Compute half-step velocity
  axpy (3, 0.5 * delta / m, f, 1, v, 1);
  printf ("fx=%le fy=%le fz=%le\n", f[0], f[1], f[2]);
  printf ("tx=%le ty=%le tz=%le\n", mol->t[0], mol->t[1], mol->t[2]);
  // Reposition center of mass
  axpy (3, delta, v, 1, com, 1);
  //printf ("Velocity norm: %le\n", nrm2 (3, v, 1));
  
  // Rotational motion
  
  // Compute angular velocity
  angular_velocity_molecule (q, qv, w);
  // Compute quaternion acceleration
  quaternion_acceleration_molecule (w, I, tp, q, qv, qa);
  // Compute half-step quaternion velocity
  axpy (4, 0.5 * delta, qa, 1, qv, 1);
  // Compute half-step quaternion as arithmetic mean
  copy_quaternion (q, qh);
  axpy (4, 0.5 * delta, qv, 1, qh, 1);
  // Compute new quaternion
  axpy (4, delta, qv, 1, q, 1);
  // Build new rotation matrix and rotate molecule
  adjust_molecule (q, mol);
}

real
velocityVerlet_molecule (real delta, pmolecule mol) {
  pquaternion q, qv, qa, qh;
  preal v, w, w_old, f, t, tp, R, I;
  real m, norm, eps, E;
  
  w = mol->w;
  I = mol->I;
  R = mol->R;
  q = mol->q;
  qv = mol->qv;
  qa = mol->qa;
  qh = mol->qh;
  v = mol->v;
  f = mol->f;
  t = mol->t;
  tp = mol->tp;
  m = mol->m;
  w_old = allocreal (3);
  
  // Translational motion
  
  // Compute new velocity
  axpy (3, 0.5 * delta / m, f, 1, v, 1);
  //printf ("||f|| = %le  ||v|| = %le\n", nrm2(3,f,1), nrm2(3,v,1));
  
  // Rotational motion
  
  // Transfer torque to reference frame
  set_zero (3, tp);
  gemv (true, 3, 3, 1.0, R, 3, t, 1, 1.0, tp, 1);
  // Use half-step angular velocity as initial approximation
  angular_velocity_molecule (qh, qv, w);
  // Save half-step quaternion velocity in intermediary storage
  copy_quaternion (qv, qh);
  // Iterate computation of angular and quaternion velocity 
  eps = 1.0;
  norm = 1.0;
  while (eps / norm > KIPS_MACH_EPS) {
    // Compute new iterate of quaternion acceleration
    quaternion_acceleration_molecule (w, I, tp, q, qv, qa);
    // Compute new iterate of quaternion velocity
    copy_quaternion (qh, qv);
    axpy (4, 0.5 * delta, qa, 1, qv, 1);
    // Save old iterate of angular velocity
    w_old[0] = w[0]; w_old[1] = w[1]; w_old[2] = w[2];
    norm = nrm2 (3, w_old, 1);
    // Compute new iterate of angular velocity
    angular_velocity_molecule (q, qv, w);
    axpy (3, -1.0, w, 1, w_old, 1);
    eps = nrm2 (3, w_old, 1);
  }
  
  // Kinetic energy
  E = 0.5 * (m * dot (false, false, 3, v, 1, v, 1) + w[0] * w[0] * I[0] 
              + w[1] * w[1] * I[1] + w[2] * w[2] * I[2]);
              
  freemem (w_old);
              
  return E;
}
