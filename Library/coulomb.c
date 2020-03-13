#include "coulomb.h"
#include <math.h>
#include <stdio.h>

/** ----------------------------------------------------------------
 *    Range splitting functions for cutoff schemes
 *  ---------------------------------------------------------------- */

real 
split_cutoff (real alpha, real R, real r) {
  (void) alpha;
  (void) R;
  
  return 1.0;
}

real
derivative_cutoff (real alpha, real R, real r) {
  (void) alpha;
  (void) R;
  
  return 0.0;
}

real 
split_SP (real alpha, real R, real r) {
  (void) alpha;
  
  return 1.0 - r/R;
}

real
derivative_SP (real alpha, real R, real r) {
  (void) alpha;
  (void) r;
  (void) R;
  
  return -1.0;
}

real 
split_SF (real alpha, real R, real r) {
  (void) alpha; 
  real z = r / R;
  
  return 1.0 - z * (2.0 - z); 
}

real 
derivative_SF (real alpha, real R, real r) {
  (void) alpha;
  real z = r / R;
  
  return z * z - 1.0;
}

real 
split_DSP (real alpha, real R, real r) {  
  return REAL_ERFC(alpha*r) - REAL_ERFC(alpha*R) * r/R;
}

real
derivative_DSP (real alpha, real R, real r) {
  return 2.0 * alpha * r / REAL_SQRT(M_PI) * REAL_EXP(-alpha*alpha*r*r) 
          - REAL_ERFC(alpha*r) - REAL_ERF(alpha*R) * r / R;
}

real 
split_DSF (real alpha, real R, real r) {
  real a, b, z;
  
  a = REAL_ERFC (alpha*R);
  b = 2.0*alpha*R/REAL_SQRT(M_PI) * REAL_EXP(-alpha*alpha*R*R);
  z = r / R;
  
  return REAL_ERFC(alpha*r) - z * (2.0*a + b - z * (a + b));
}

real 
derivative_DSF (real alpha, real R, real r) {
  real a, b, z;
  
  a = REAL_ERFC (alpha*R);
  b = 2.0 * alpha * R / REAL_SQRT(M_PI) * REAL_EXP(-alpha*alpha*R*R);
  z = r / R;
  return (a + b) * z * z - REAL_ERFC (alpha*r) 
          - 2.0 * alpha * r / REAL_SQRT(M_PI) * REAL_EXP(-alpha*alpha*r*r);
}

real 
split_SP2 (real alpha, real R, real r) {
  real z = r / R;
  (void) alpha;
  
  return 1.0 - z * (2.0 - z*z * (2.0 - z)); 
}

real
derivative_SP2 (real alpha, real R, real r) {
  real z = r / R;
  (void) alpha;
  
  return -1.0 + z * z * z * (4.0 - 3.0 * z);
}

real 
split_SP3 (real alpha, real R, real r) {
  real z = r / R;
  real y = z * z;
  (void) alpha;
  
  return 1.0 - z * (1.75 - y * y * (5.25 - z * (7.0 - 2.5*z)));
}

real 
derivative_SP3 (real alpha, real R, real r) {
  real z = r / R;
  real y = z * z;
  (void) alpha;
  
  return -1.0 + y * y * z * (21.0 - z * (35.0 - 15.0 * z));
}


/** ----------------------------------------------------------------
 *    Geometric data
 *  ---------------------------------------------------------------- */

real
mindiam_coulomb (real eta, real rmol) {
  return 2.0 * eta * rmol;
}


/** ----------------------------------------------------------------
 *    Parameters
 *  ---------------------------------------------------------------- */
 
pcoulomb
setParametersPotential_coulomb (split s, real alpha, real R, pcspatialgeometry sg) {
  pcoulomb par = (pcoulomb) allocmem (sizeof (coulomb));
  par->s = s;
  par->alpha = alpha;
  par->R = R;
  par->sg = sg;
  
  return par;
}

pcoulomb
setParametersGradient_coulomb (split s, split sd, real alpha, real R, 
                               pcspatialgeometry sg) {
  pcoulomb par = (pcoulomb) allocmem (sizeof (coulomb));
  par->s = s;
  par->sd = sd;
  par->alpha = alpha;
  par->R = R;
  par->sg = sg;
  
  return par;
}

void
delParameters_coulomb (pcoulomb c) {
  freemem (c);
}

/** ----------------------------------------------------------------
 *    Coulomb Potential
 *  ---------------------------------------------------------------- */
 
/** Evaluate the pairwise potential between all periodic copies of x and y within 
 *  the cutoff distance R, recursively by dimension. */
static real 
sum_potential_cutoff (uint dim, uint ldim, uint *j, uint *m, pcreal d, int *j0, 
                      bool exclude, split s, real alpha, real R, pcreal x, pcreal y){
  int i;
  real a, b, c, R2, r, V;
  V = 0.0;
  bool exclude_new;
  
  if (dim > 1) {
    dim --;
    for (i=0; i<m[dim]; i++){
      j[dim] = i;
      if (i == j0[dim]) {
        exclude_new = exclude;
      }
      else {
        exclude_new = false;
      }
      V += sum_potential_cutoff (dim, ldim, j, m, d, j0, exclude_new, s, alpha, R, x, y);
    }
  }
  else {
    assert (dim == 1);
    R2 = R*R;
    c = 0.0;
    for (i=ldim-1; i>=0; i--){
      a = y[i] + j[i] * d[i] - x[i];
      b = a * a;
      c += b;
    }
    for (i=0; i<m[0]; i++) { 
      if (i != j0[0] || exclude == false) { 
        if (c < R2) {
          r = REAL_SQRT (c);
          V += s (alpha, R, r) / r;
        }
      }
      c -= b;
      a += d[0];
      b = a * a;
      c += b;
    }
  }
  return V;
}

static real 
shortrangePotential_coulomb (pcreal x, pcreal y, bool exclude, pcoulomb par) {
  uint i, dim, *j, *m;
  int jmin, jmax, *j0;
  real V, R, alpha;
  preal d, z;
  pcspatialgeometry sg = par->sg;
  split s = par->s;
  
  alpha = par->alpha;
  dim = sg->dim;
  R = par->R;
  z = allocreal (dim);
  d = allocreal (dim);
  j = allocuint (dim);
  j0 = (int *) allocmem (sizeof (int) * dim);
  m = allocuint (dim);
  
  // Determine all periodic copies of y within the given radius around x.
  for (i=0; i<dim; i++) {
    d[i] = sg->bmax[i] - sg->bmin[i];
    jmin = ceil ((x[i] - y[i] - R) / d[i]);
    jmax = floor ((x[i] - y[i] + R) / d[i]);
    m[i] = (jmax < jmin ? 0 : jmax-jmin+1);
    z[i] = y[i] + jmin * d[i];
    j[i] = 0;
    j0[i] = -jmin;
  }
  
  V = sum_potential_cutoff (dim, dim, j, m, d, j0, exclude, s, alpha, R, x, z);
  
  freemem (j);
  freemem (m);
  freemem (d);
  freemem (z);
  freemem (j0);
  
  return V;
}


real 
kernel_coulomb (pcreal x, pcreal y, uint xmol, uint ymol, void *data) {
  real V;
  
  /** If the molecule numbers of x and y are the same and they are not 
   *  interpolation points, then their interaction will be neglected. 
   *  Interactions with periodic copies will nevertheless be considered. */
  if (xmol == ymol && xmol > 0) {
    V = shortrangePotential_coulomb (x, y, true, (pcoulomb) data);
  }
  else {
    V = shortrangePotential_coulomb (x, y, false, (pcoulomb) data);
  }
  return V;
}

real 
selfterm_coulomb (split sp, real alpha, pcavector q, real R) {
  real V;
  
  V = -REAL_SQR (norm2_avector (q))/R;
  
  if (sp == &split_cutoff) {
    V *= 0.0;
  }
  else if (sp == &split_SP) {
    V *= 0.5;
  }
  else if (sp == &split_DSP) {
    V *= (0.5 * REAL_ERFC (alpha * R) + R * alpha/REAL_SQRT(M_PI));
  }
  else if (sp == &split_DSF) {
    V *= (REAL_ERFC (alpha * R) + R * alpha/REAL_SQRT(M_PI) * 
          (1.0 + REAL_EXP (- alpha * alpha * R * R)));
  }
  else if (sp == &split_SP3) {
    V *= 7.0/8.0;
  }
  
  return V;
}

real
energy_coulomb (pcavector q, pckernelmatrix kc, pspatialgeometry sg, ph2matrix Vc) {
  real E;
  real C;
  uint n = kc->points;
  pavector y = new_avector (n);
  pcoulomb pc = (pcoulomb) kc->data;
  
  clear_avector (y);
  assert (n == q->size);
  
  C = Q_E * Q_E / (4.0 * M_PI * EPS_0);
  initPoints_spatialgeometry (n, kc->x, sg);
  update_h2matrix_kernelmatrix (false, kc, Vc);
  mvm_h2matrix (1.0, false, Vc, q, y);
  E = C * (0.5 * dotprod_avector (q, y) + selfterm_coulomb (pc->s, pc->alpha, q, pc->R));
  
  del_avector (y);
  
  return E;
}

/** ----------------------------------------------------------------
 *    Coulomb Forces
 *  ---------------------------------------------------------------- */
 
/** Evaluate the gradient of the pairwise potential between all periodic copies of x 
 *  and y within the cutoff distance R, recursively by dimension. */
static void
sum_gradient_cutoff (uint dim, uint ldim, uint *j, uint *m, pcreal d, int *j0, 
                     bool exclude, split s, split sd, real alpha, real R, pcreal x, 
                     pcreal y, pfield f){
  int i;
  uint l;
  real a, b, r2, R2, r, z[ldim];
  bool exclude_new;
  
  if (dim > 1) {
    dim --;
    for (i=0; i<m[dim]; i++){
      j[dim] = i;
      if (i == j0[dim]) {
        exclude_new = exclude;
      }
      else {
        exclude_new = false;
      }
      sum_gradient_cutoff (dim, ldim, j, m, d, j0, exclude_new, s, sd, alpha, R, x, y, f);
    }
  }
  else {
    assert (dim == 1);
    R2 = R*R;
    r2 = 0.0;
    for (i=ldim-1; i>=0; i--){
      z[i] = y[i] + j[i] * d[i] - x[i];
      b = z[i] * z[i];
      r2 += b;
    }
    for (i=0; i<m[0]; i++) {
      if (i != j0[0] || exclude == false) {
        if (r2 < R2) {
          r = REAL_SQRT (r2);
          a = (s (alpha, r, R) - sd (alpha, r, R) * r) / (r * r2);
          for (l=0; l<ldim; l++) {
            f[l] += z[l] * a;
          }
        }
      }  
      r2 -= b;
      z[0] += d[0];
      b = z[0] * z[0];
      r2 += b;
    }
  }
}

static void
shortrangeGradient_coulomb (pcreal x, pcreal y, bool exclude, pcoulomb par, pfield f) {
  uint i, dim, *j, *m;
  int jmin, jmax, *j0;
  preal d, z;
  real alpha, R;
  pcspatialgeometry sg = par->sg;
  split s, sd;
  
  alpha = par->alpha;
  R = par->R;
  s = par->s;
  sd = par->sd;
  dim = sg->dim;
  z = allocreal (dim);
  d = allocreal (dim);
  j = allocuint (dim);
  m = allocuint (dim);
  j0 = (int *) allocmem (sizeof (int) * dim);
  
  // Determine all periodic copies of y within the given radius around x.
  for (i=0; i<dim; i++) {
    d[i] = sg->bmax[i] - sg->bmin[i];
    jmin = ceil ((x[i] - y[i] - R) / d[i]);
    jmax = floor ((x[i] - y[i] + R) / d[i]);
    m[i] = (jmax < jmin ? 0 : jmax-jmin+1);
    z[i] = y[i] + jmin * d[i];
    j[i] = 0;
    f[i] = 0.0;
    j0[i] = jmin;
  }
  
  sum_gradient_cutoff (dim, dim, j, m, d, j0, exclude, s, sd, alpha, R, x, z, f);
  
  freemem (j);
  freemem (m);
  freemem (d);
  freemem (z);
  freemem (j0);
}

void 
gradient_coulomb (pcreal x, pcreal y, uint xmol, uint ymol, void *data, pfield f) {
  /** If the molecule numbers of x and y are the same and they are not 
   *  interpolation points, then their interaction will be neglected. 
   *  Interactions with periodic copies will nevertheless be considered. */
  if (xmol == ymol && xmol > 0) {
    shortrangeGradient_coulomb (x, y, true, (pcoulomb) data, f);
  }
  else {
    shortrangeGradient_coulomb (x, y, false, (pcoulomb) data, f);
  }
}

void
fullForce_coulomb (pcavector q, pckernelmatrix kc, pspatialgeometry sg, 
                   ph2matrix Fc, pavector f) {
  uint dim = kc->dim;
  uint points = kc->points;
  uint i, j, off;
  real C;
  
  C = Q_E * Q_E / (4.0 * M_PI * EPS_0);
  
  assert (q->size == points);
  assert (f->size == points * dim);
  clear_avector (f);
  initPoints_spatialgeometry (points, kc->x, sg);
  update_h2matrix_kernelmatrix (true, kc, Fc);
  mvm_h2matrix (1.0, false, Fc, q, f);
  for (i=0; i<points; i++) {
    off = i * dim;
    for (j=0; j<dim; j++) {
      f->v[j+off] *= -q->v[i];
    }
  }
  scale_avector (C, f);
}
