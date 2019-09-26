#include "coulomb.h"
#include <math.h>
#include <omp.h>
#include <stdio.h>

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
  return 2.0 * alpha * r / REAL_SQRT(M_PI) * REAL_EXP(-alpha*alpha*r*r) - REAL_ERFC(alpha*r) 
          - REAL_ERF(alpha*R) * r / R;
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

static real 
sum_potential_cutoff (uint dim, uint ldim, uint *j, uint *m, pcreal d, split s, 
                      real alpha, real R, pcreal x, pcreal y){
  int i;
  real a, b, c, R2, r, V;
  V = b = 0.0;
  
  assert (R > KIPS_ALMOST_ZERO);
  
  if (dim > 1) {
    dim --;
    for (i=0; i<m[dim]; i++){
      j[dim] = i;
      V += sum_potential_cutoff (dim, ldim, j, m, d, s, alpha, R, x, y);
    }
  }
  else {
    assert (dim == 1);
    R2 = R*R;
    c = 0.0;
    for (i=ldim-1; i>=0; i--){
      a = x[i] + j[i] * d[i] - y[i];
      b = a * a;
      c += b;
    }
    for (i=0; i<m[0]; i++) {
      if (c > KIPS_MACH_EPS && c < R2) {
        r = REAL_SQRT (c);
        V += s (alpha, R, r) / r;
      }
      c -= b;
      a += d[0];
      b = a * a;
      c += b;
    }
  }
  return V;
}

real 
shortrangePotential_coulomb (split s,real alpha, real R, pspatialgeometry sg, 
                             pcreal x, pcreal y) {
  uint i, dim, *j, *m;
  int jmin, jmax;
  real V;
  preal d, z;
  
  dim = sg->dim;
  z = allocreal (dim);
  d = allocreal (dim);
  j = allocuint (dim);
  m = allocuint (dim);
  
  for (i=0; i<dim; i++) {
    d[i] = sg->bmax[i] - sg->bmin[i];
    jmin = ceil ((y[i] - x[i] - R)/d[i]);
    jmax = floor ((y[i] - x[i] + R)/d[i]);
    m[i] = (jmax < jmin ? 0 : jmax-jmin+1);
    z[i] = x[i] + jmin * d[i];
    j[i] = 0;
  }
  
  V =sum_potential_cutoff (dim, dim, j, m, d, s, alpha, R, z, y);
  
  freemem (j);
  freemem (m);
  freemem (d);
  freemem (z);
  
  return V;
}

pparameters
setparametersPotential_coulomb (split s, real alpha, real R, pspatialgeometry sg) {
  pparameters p;
  
  p = (pparameters) allocmem (sizeof (parameters));
  p->s = s;
  p->alpha = alpha;
  p->R = R;
  p->sg = sg;
  
  return p;
}

real 
kernel_coulomb (pcreal x, pcreal y, void *data) {
  pparameters p = (pparameters) data;
  
  return shortrangePotential_coulomb (p->s, p->alpha, p->R, p->sg, x, y);
}

static void
sum_gradient_cutoff (uint dim, uint ldim, uint *j, uint *m, pcreal d, split s, split sd, 
                     real alpha, real R, pcreal x, pcreal y, pfield f){
  int i;
  uint l;
  real a, b, c, R2, r, z[dim];
  b = 0.0;
  
  assert (R > KIPS_ALMOST_ZERO);
  
  if (dim > 1) {
    dim --;
    for (i=0; i<m[dim]; i++){
      j[dim] = i;
      sum_gradient_cutoff (dim, ldim, j, m, d, s, sd, alpha, R, x, y, f);
    }
  }
  else {
    assert (dim == 1);
    R2 = R*R;
    c = 0.0;
    for (i=ldim-1; i>=0; i--){
      z[i] = x[i] + j[i] * d[i] - y[i];
      b = z[i] * z[i];
      c += b;
    }
    for (i=0; i<m[0]; i++) {
      if (c > KIPS_MACH_EPS && c < R2) {
        r = REAL_SQRT (c);
        a = (sd (alpha, r, R) * r - s (alpha, r, R)) / (r * c);
        for (l=0; l<dim; l++) {
          f[l] += z[l] * a;
        }
      }
      c -= b;
      z[0] += d[0];
      b = z[0] * z[0];
      c += b;
    }
  }
}

void
shortrangeGradient_coulomb (split s, split sd, real alpha, real R, pspatialgeometry sg, 
                            pcreal x, pcreal y, pfield f) {
  uint i, dim, *j, *m;
  int jmin, jmax;
  preal d, z;
  
  dim = sg->dim;
  z = allocreal (dim);
  d = allocreal (dim);
  j = allocuint (dim);
  m = allocuint (dim);
  
  for (i=0; i<dim; i++) {
    d[i] = sg->bmax[i] - sg->bmin[i];
    jmin = ceil ((y[i] - x[i] - R)/d[i]);
    jmax = floor ((y[i] - x[i] + R)/d[i]);
    m[i] = (jmax < jmin ? 0 : jmax-jmin+1);
    z[i] = x[i] + jmin * d[i];
    j[i] = 0;
    f[i] = 0.0;
  }
  
  sum_gradient_cutoff (dim, dim, j, m, d, s, sd, alpha, R, z, y, f);
  
  freemem (j);
  freemem (m);
  freemem (d);
  freemem (z);
}

pparameters
setparametersGradient_coulomb (split s, split sd, real alpha, real R, pspatialgeometry sg) {
  pparameters p;
  
  p = (pparameters) allocmem (sizeof (parameters));
  p->s = s;
  p->sd = sd;
  p->alpha = alpha;
  p->R = R;
  p->sg = sg;
  
  return p;
}

void 
gradient_coulomb (pcreal x, pcreal y, void *data, pfield f) {
  pparameters p = (pparameters) data;
  
  shortrangeGradient_coulomb (p->s, p->sd, p->alpha, p->R, p->sg, x, y, f);
}

real 
selfPotential_coulomb (split sp, real alpha, pavector q, real R) {
  real V;
  
  V = -REAL_SQR (norm2_avector (q))/R;
  
  if (sp == &split_SP) {
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
