#include "lj.h"
#include <math.h>
#include <stdio.h>

/** ----------------------------------------------------------------
 *    Parameters
 *  ---------------------------------------------------------------- */

plj
setParameters_lj (real R, real sig, real eps, pcspatialgeometry sg) {
  real sig2, sig6, sig12, R2, R6, R12;

  plj par = (plj) allocmem (sizeof (lj));
  R2 = R * R;
  R6 = R2 * R2 * R2;
  R12 = R6 * R6;
  sig2 = sig * sig;
  sig6 = sig2 * sig2 * sig2;
  sig12 = sig6 * sig6;
  par->R = R;
  par->eps = eps;
  par->sig6 = sig6;
  par->sig12 = sig12;
  par->sg = sg;
  
  // Shift parameters.
  par->a = (6.0 * sig12 / R12 - 3.0 * sig6 / R6) / R2;
  par->b = 4.0 * sig6 / R6 - 7.0 * sig12 / R12;
  
  return par;
}

plj
setParametersUnshifted_lj (real R, real sig, real eps, pcspatialgeometry sg) {
  plj par = setParameters_lj (R, sig, eps, sg);
  
  // Do not use shift.
  par->a = 0.0;
  par->b = 0.0;
  
  return par;
}

void 
delParameters_lj (plj lj) {
  freemem (lj);
}

/** ----------------------------------------------------------------
 *    Lennard-Jones potential
 *  ---------------------------------------------------------------- */

/** Evaluate the pairwise potential between all periodic copies of x and y 
 *  within the cutoff distance R, recursively by dimension. */
static real 
sum_potential_cutoff (uint dim, uint ldim, uint *j, uint *m, pcreal d, 
                      int *j0, bool exclude, real R, real sig12, real sig6, 
                      real alpha, real beta, pcreal x, pcreal y){
  int i;
  real a, b, r2, R2, V, r6, r12;
  bool exclude_new;
  V = 0.0;
  
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
      V += sum_potential_cutoff (dim, ldim, j, m, d, j0, exclude_new, R, sig12, 
                                 sig6, alpha, beta, x, y);
    }
  }
  else {
    assert (dim == 1);
    R2 = R * R;
    r2 = 0.0;
    for (i=ldim-1; i>=0; i--){
      a = y[i] + j[i] * d[i] - x[i];
      b = a * a;
      r2 += b;
    }
    for (i=0; i<m[0]; i++) {
      if (i != j0[0] || exclude == false) {
        if (r2 < R2) {
          r6 = r2 * r2 * r2;
          r12 = r6 * r6;
          V += sig12 / r12 - sig6 / r6 + alpha * r2 + beta;
        }
      }
      r2 -= b;
      a += d[0];
      b = a * a;
      r2 += b;
    }
  }
  return V;
}

static real
shortrangePotential_lj (pcreal x, pcreal y, bool exclude, plj par) {
  uint i, dim, *j, *m;
  int jmin, jmax, *j0;
  real V, R;
  preal d, z;
  pcspatialgeometry sg = par->sg;
  R = par->R;
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
    j0[i] = jmin;
  }
  
  V = 4.0 * par->eps;
  V *= sum_potential_cutoff (dim, dim, j, m, d, j0, exclude, R, par->sig12, 
                              par->sig6, par->a, par->b, x, z);
  
  freemem (j);
  freemem (j0);
  freemem (m);
  freemem (d);
  freemem (z);
  
  return V;
}

real
kernel_lj (pcreal x, pcreal y, uint xmol, uint ymol, void *data) {
  real V;
  
  /** If the molecule numbers of x and y are the same and they are not 
   *  interpolation points, then their interaction will be neglected. 
   *  Interactions with periodic copies will nevertheless be considered. */
  if (xmol == ymol && xmol > 0) {
    V = shortrangePotential_lj (x, y, true, (plj) data);
  }
  else {
    V = shortrangePotential_lj (x, y, false, (plj) data);
  }
  
  return V;
}

real
energy_lj (pckernelmatrix klj, pspatialgeometry sg, ph2matrix Vlj) {
  real E;
  pavector x, z;
  uint n = klj->points;
  
  x = new_avector (n);
  z = new_avector (n);
  fill_avector (1.0, x);
  clear_avector (z);
  
  initPoints_spatialgeometry (n, klj->x, sg);
  update_h2matrix_kernelmatrix (false, klj, Vlj);
  mvm_h2matrix (1.0, false, Vlj, x, z);
  E = 0.5 * dotprod_avector (x, z);
  
  del_avector (x);
  del_avector (z);
  
  return E;
}


/** ----------------------------------------------------------------
 *    Lennard-Jones forces
 *  ---------------------------------------------------------------- */

/** Evaluate the gradient of the pairwise potential between all periodic copies 
 *  of x and y within the cutoff distance R, recursively by dimension. */
static void
sum_gradient_cutoff (uint dim, uint ldim, uint *j, uint *m, pcreal d, int *j0, 
                     bool exclude, real R, real sig12, real sig6, real alpha, 
                     pcreal x, pcreal y, pfield f){
  int i;
  uint l;
  real a, b, r2, R2, r4, r8, r16, z[ldim];
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
      sum_gradient_cutoff (dim, ldim, j, m, d, j0, exclude_new, R, sig12, sig6, alpha, 
                           x, y, f);
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
        //printf ("r^2=%le\n", r2);
        if (r2 < R2) {
          r4 = r2 * r2;
          r8 = r4 * r4;
          r16 = r8 * r8;
          a = 12.0 * sig12 * r2 / r16 - 6.0 * sig6 / r8 + 2.0 * alpha;
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
shortrangeGradient_lj (pcreal x, pcreal y, bool exclude, plj par, pfield f) {
  uint i, dim, *j, *m;
  int jmin, jmax, *j0;
  preal d, z;
  real R, eps;
  pcspatialgeometry sg = par->sg;
  
  R = par->R;
  eps = par->eps;
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
    //printf ("%d  %d  %u  %u\n", jmin, jmax, m[i], exclude);
  }
  
  sum_gradient_cutoff (dim, dim, j, m, d, j0, exclude, R, par->sig12, par->sig6, 
                        par->a, x, z, f);
  for (i=0; i<dim; i++) {
    f[i] *= 4.0 * eps;
  }
  
  freemem (j);
  freemem (m);
  freemem (d);
  freemem (z);
  freemem (j0);
}

void
gradient_lj (pcreal x, pcreal y, uint xmol, uint ymol, void *data, pfield f) {
  if (xmol == ymol && xmol > 0) {
    shortrangeGradient_lj (x, y, true, (plj) data, f);
  }
  else {
    shortrangeGradient_lj (x, y, false, (plj) data, f);
  }
  //printf ("%le\n", sqrt((x[0]-y[0])*(x[0]-y[0])+(x[1]-y[1])*(x[1]-y[1])+(x[2]-y[2])*(x[2]-y[2])));
}

void
fullForce_lj (pckernelmatrix klj, pspatialgeometry sg, ph2matrix Flj, pavector f) {
  pavector q;
  uint dim = klj->dim;
  uint points = klj->points;
  
  assert (f->size == dim * points);
  q = new_avector (points);
  fill_avector (-1.0, q);
  clear_avector (f);
  
  initPoints_spatialgeometry (points, klj->x, sg);
  update_h2matrix_kernelmatrix (true, klj, Flj);
  mvm_h2matrix (1.0, false, Flj, q, f);
  del_avector (q);
}
