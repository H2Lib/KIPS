#include "coulomb.h"
#include <math.h>
#include <omp.h>
#include <stdio.h>

real 
split_SP (real alpha, real R, real r) {
  (void) alpha;
  
  assert (R > KIPS_ALMOST_ZERO);
  
  return 1.0 - r/R;
}

real 
split_SF (real alpha, real R, real r) {
  real z;;
  (void) alpha;
  
  assert (R > KIPS_ALMOST_ZERO);
  z = r/R;
  
  return 1.0 - z * (2.0 - z); 
}

real 
split_DSP (real alpha, real R, real r) {  
  assert (R > KIPS_ALMOST_ZERO);

  return erfc (alpha*r) - erfc (alpha*R) * r/R;
}

real 
split_DSF (real alpha, real R, real r) {
  real a, b, z;
  
  assert (R > KIPS_ALMOST_ZERO);
  a = erfc (alpha*R);
  b = 2.0*alpha/sqrt(M_PI) * exp (-REAL_SQR(alpha)*REAL_SQR(R))/R;
  z = r/R;
  
  return erfc (alpha*r) - z * (2.0*a + b - z * (a + b));
}

real 
split_SP2 (real alpha, real R, real r) {
  real z;
  (void) alpha;
  
  assert (R > KIPS_ALMOST_ZERO);
  z = r/R;
  
  return 1.0 - z * (2.0 - z*z * (2.0 - z)); 
}

real 
split_SP3 (real alpha, real R, real r) {
  real z;
  (void) alpha;
  
  assert (R > KIPS_ALMOST_ZERO);
  z = r/R;
  
  return 1.0 -  z * (1.75 - z*z*z*z * (5.25 - z * (7.0 - 2.5*z)));
}

static real 
sumPotential_coulomb (uint dim, uint ldim, uint *j, uint *m, pcreal d, split s, real alpha, real R, pcreal x, pcreal y){
  int i;
  real a, b, c, R2, r, V;
  V = b = 0.0;
  
  if (dim > 1) {
    for (i=0; i<m[dim-1]; i++){
      j[dim-1] = i;
      V += sumPotential_coulomb (dim-1, ldim, j, m, d, s, alpha, R, x, y);
    }
  }
  else {
    assert (dim == 1);
    R2 = R*R;
    c = 0.0;
    for (i=ldim-1; i>=0; i--){
      a = x[i]+j[i]*d[i]-y[i];
      b = a*a;
      c += b;
    }
    for (i=0; i<m[0]; i++) {
      if (c < R2 && c > H2_MACH_EPS) {
        r = sqrt(c);
        V += s(alpha,R,r)/r;
      }
      c -= b;
      a += d[0];
      b = a*a;
      c += b;
    }
  }
  return V;
}

real 
shortrangePotential_coulomb (split s, real alpha, real R, pspatialgeometry sg, pcreal x, pcreal y) {
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
  
  V =sumPotential_coulomb (dim, dim, j, m, d, s, alpha, R, z, y);
  
  freemem (j);
  freemem (m);
  freemem (d);
  freemem (z);
  
  return V;
}

pparameters
setparameters_coulomb (split s, real alpha, real R, pspatialgeometry sg) {
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

real 
selfPotential_coulomb (split sp, real alpha, pavector q, real R) {
  real V;
  
  V = -REAL_SQR (norm2_avector (q))/R;
  
  if (sp == &split_SP) {
    V *= 0.5;
  }
  else if (sp == &split_DSP) {
    V *= (0.5 * erfc (alpha * R) + R * alpha/REAL_SQRT(M_PI));
  }
  else if (sp == &split_DSF) {
    V *= (erfc (alpha * R) + R * alpha/REAL_SQRT(M_PI) * (1.0 + exp (- alpha * alpha * R * R)));
  }
  else if (sp == &split_SP3) {
    V *= 7.0/8.0;
  }
  
  return V;
}
