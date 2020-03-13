#include "tip4p.h"
#include <stdio.h>

/** ---------------------------------------------------------------------
 *    Constructors and destructors
 *  --------------------------------------------------------------------- */

ptip4p 
new_tip4p (uint n) {
  uint i, j;
  preal *xc, *xn;
  pmolecule *mol;
  pfield q;
  ptip4p t = (ptip4p) allocmem (sizeof (tip4p));
  
  t->n = n;
  t->mol = mol = (pmolecule *) allocmem (n * sizeof (pmolecule));
  t->xc = xc = (preal *) allocmem (3 * n * sizeof (preal));
  t->xn = xn = (preal *) allocmem (n * sizeof (preal));
  
  xc[0] = allocreal (3 * n * 3);
  xn[0] = allocreal (n * 3);
  
  t->xcmol = allocuint (3 * n);
  t->xnmol = allocuint (n);
  t->q = new_avector (3 * n);
  
  q = t->q->v;
  for (i=0; i<n; i++) {
    j = 3 * i;
    mol[i] = new_molecule (3);
    xn[i] = xn[0] + j;
    xc[j] = xc[0] + 3 * j;
    xc[j+1] = xc[j] + 3;
    xc[j+2] = xc[j] + 6;
    q[j] = QM_TIP4P;
    q[j+1] = q[j+2] = QH_TIP4P;
  }
  
  return t;
}

void
del_tip4p (ptip4p t) {
  for (uint i=0; i<t->n; i++) {
    del_molecule (t->mol[i]);
  }
  freemem (t->xc[0]);
  freemem (t->xc);
  freemem (t->xn[0]);
  freemem (t->xn);
  freemem (t->xcmol);
  freemem (t->xnmol);
  del_avector (t->q);
  freemem (t->mol);
  freemem (t);
}
 
/** ---------------------------------------------------------------------
 *    In-/Output
 *  --------------------------------------------------------------------- */

static void 
calc_virtual_charge (pcreal h1, pcreal h2, pcreal o, preal m) {
  real xm, ym, zm, rm;
  real rOM = ROM_TIP4P;
  
  xm = 0.5 * (h1[0] + h2[0]) - o[0];
  ym = 0.5 * (h1[1] + h2[1]) - o[1];
  zm = 0.5 * (h1[2] + h2[2]) - o[2];
  rm = rOM / REAL_SQRT (xm * xm + ym * ym + zm * zm);
  m[0] = o[0] + xm * rm;
  m[1] = o[1] + ym * rm;
  m[2] = o[2] + zm * rm;
}

ptip4p 
inputfile_tip4p (const char *file) {
  ptip4p t;
  FILE *stream;
  uint a, b, c, i, j, k, l, n, N, *xnmol, *xcmol;
  bool *sw;
  char buf[100];
  preal *xn, *xc, r, com;
  pmolecule *mol;
  real q, x, y, z;
  real M = 2 * M_H + M_O;
  stream = fopen (file, "r");
  
  if (stream == NULL) {
    printf ("Error: File not found.\n");
    exit (0);
  }
  fscanf (stream, "%u", &n);
  fgets (buf, 100, stream);
  
  assert (n%3 == 0);
  N = n / 3;
  t = new_tip4p (N);
  xn = t->xn;
  xc = t->xc;
  xnmol = t->xnmol;
  xcmol = t->xcmol;
  mol = t->mol;
  sw = allocmem (N * sizeof (bool));
  for (i=0; i<N; i++) {
    sw[i] = false;
  }
  for (i=0; i<n; i++) {
    #ifdef USE_FLOAT
      fscanf (stream, "%u %u %u %e %e %e %e %u %u %u\n", &j, &k, &l, &q, &x, &y, &z, &a, &b, &c);
    #else
      fscanf (stream, "%u %u %u %le %le %le %le %u %u %u\n", &j, &k, &l, &q, &x, &y, &z, &a, &b, &c);
    #endif
    if (l == 1) {
      xn[k-1][0] = x; xn[k-1][1] = y; xn[k-1][2] = z;
      xnmol[k-1] = k;
    }
    else if (l == 2) {
      j = 3 * (k - 1) + 1;
      if (sw[k-1] == false){
        sw[k-1] = true;
      }
      else {
        j += 1;
      }
      xc[j][0] = x; xc[j][1] = y; xc[j][2] = z;
      xcmol[j] = k;
    }
    else {
      printf ("Error: Unknown atom type %u in line %u.\n", l, i+3);
      exit (0);
    }
  }
  for (i=0; i<N; i++) {
    j = 3 * i;
    // Calculate virtual charge location.
    calc_virtual_charge (xc[j+1], xc[j+2], xn[i], xc[j]);
    xcmol[j] = i + 1;
   
    // Calculate center of mass.
    com = mol[i]->com;
    com[0] = (M_H * (xc[j+1][0] + xc[j+2][0]) + M_O * xn[i][0]) / M;
    com[1] = (M_H * (xc[j+1][1] + xc[j+2][1]) + M_O * xn[i][1]) / M;
    com[2] = (M_H * (xc[j+1][2] + xc[j+2][2]) + M_O * xn[i][2]) / M;
    mol[i]->m = M;
    
    // Calculate relative atom positions.
    r = mol[i]->r[0];
    r[0] = xn[i][0] - com[0];
    r[1] = xn[i][1] - com[1];
    r[2] = xn[i][2] - com[2];
    r = mol[i]->r[1];
    r[0] = xc[j+1][0] - com[0];
    r[1] = xc[j+1][1] - com[1];
    r[2] = xc[j+1][2] - com[2];
    r = mol[i]->r[2];
    r[0] = xc[j+2][0] - com[0];
    r[1] = xc[j+2][1] - com[1];
    r[2] = xc[j+2][2] - com[2];
  }
  
  fclose (stream);
  freemem (sw);
  
  return t;
}

/*
ptip4p 
inputfile_tip4p (const char *file) {
  ptip4p t;
  FILE *stream;
  uint i, j, n, *xnmol, *xcmol;
  char s[3];
  char buf[100];
  preal *xn, *xc;
  pfield r, com;
  pmolecule *mol;
  real xh1, xh2, yh1, yh2, zh1, zh2, xo, yo, zo;
  real M = 2 * M_H + M_O;
  stream = fopen (file, "r");
  
  if (stream == NULL) {
    printf ("Error: File not found.\n");
    exit (0);
  }
  fscanf (stream, "%u\n", &n);
  assert (n%3 == 0);
  n /= 3;
  fgets (buf, 100, stream);
  t = new_tip4p (n);
  xn = t->xn;
  xc = t->xc;
  xnmol = t->xnmol;
  xcmol = t->xcmol;
  mol = t->mol;
  for (i=0; i<n; i++) {
    j = 3*i;
    #ifdef USE_FLOAT
    {
      fscanf (stream, "%s %f %f %f\n", s, &xo, &yo, &zo);
    }
    #else 
    {
      fscanf (stream, "%s %lf %lf %lf\n", s, &xo, &yo, &zo);
    }
    #endif
    if (s[0] != 'O' || s[1] != '\0') {
      printf ("Error: Expected oxygen atom in line %u\n", i+2);
      exit (0);
    }
    xn[i][0] = 1e-10 * xo;
    xn[i][1] = 1e-10 * yo;
    xn[i][2] = 1e-10 * zo;
    xnmol[i] = i+1;
    
    #ifdef USE_FLOAT
    {
      fscanf (stream, "%s %f %f %f\n", s, &xh1, &yh1, &zh1);
    }
    #else 
    {
      fscanf (stream, "%s %lf %lf %lf\n", s, &xh1, &yh1, &zh1);
    }
    #endif
    if (s[0] != 'H' || s[1] != '\0') {
      printf ("Error: Expected hydrogen atom in line %u\n", i+3);
      exit (0);
    }
    xc[j+1][0] = 1e-10 * xh1;
    xc[j+1][1] = 1e-10 * yh1;
    xc[j+1][2] = 1e-10 * zh1;
    xcmol[j+1] = i+1;
    
    #ifdef USE_FLOAT
    {
      fscanf (stream, "%s %f %f %f\n", s, &xh2, &yh2, &zh2);
    }
    #else
    {
      fscanf (stream, "%s %lf %lf %lf\n", s, &xh2, &yh2, &zh2);
    }
    #endif
    if (s[0] != 'H' || s[1] != '\0') {
      printf ("Error: Expected hydrogen atom in line %u\n", i+4);
      exit (0);
    }
    xc[j+2][0] = 1e-10 * xh2;
    xc[j+2][1] = 1e-10 * yh2;
    xc[j+2][2] = 1e-10 * zh2;
    xcmol[j+2] = i+1;
    
    // Calculate virtual charge location.
    calc_virtual_charge (xc[j+1], xc[j+2], xn[i], xc[j]);
    xcmol[j] = i+1;
    
    // Calculate center of mass.
    com = mol[i]->com->v;
    com[0] = (M_H * (xc[j+1][0] + xc[j+2][0]) + M_O * xn[i][0]) / M;
    com[1] = (M_H * (xc[j+1][1] + xc[j+2][1]) + M_O * xn[i][1]) / M;
    com[2] = (M_H * (xc[j+1][2] + xc[j+2][2]) + M_O * xn[i][2]) / M;
    mol[i]->m = M;
    
    // Calculate relative atom positions.
    r = mol[i]->r[0]->v;
    r[0] = xn[i][0] - com[0];
    r[1] = xn[i][1] - com[1];
    r[2] = xn[i][2] - com[2];
    r = mol[i]->r[1]->v;
    r[0] = xc[j+1][0] - com[0];
    r[1] = xc[j+1][1] - com[1];
    r[2] = xc[j+1][2] - com[2];
    r = mol[i]->r[2]->v;
    r[0] = xc[j+2][0] - com[0];
    r[1] = xc[j+2][1] - com[1];
    r[2] = xc[j+2][2] - com[2];
  }
  
  fclose (stream);
  
  return t;
}
*/

/** ---------------------------------------------------------------------
 *    Energy and force calculation
 *  --------------------------------------------------------------------- */

real
energy_tip4p (pspatialgeometry sg, pkernelmatrix kc, pkernelmatrix klj, 
                 ph2matrix Vc, ph2matrix Vlj, ptip4p t) {
  real E;
  
  update_kernelmatrix (3 * t->n, t->xc, t->xcmol, kc);
  update_kernelmatrix (t->n, t->xn, t->xnmol, klj);
  E = energy_coulomb (t->q, kc, sg, Vc) + energy_lj (klj, sg, Vlj);
  
  return E;
}

static void 
set_zero (uint n, preal x) {
  for (uint i=0; i<n; i++) 
    x[i] = 0.0;
}

real
calcForce_tip4p (pavector fc, pavector flj, ptip4p t) {
  preal f;
  pfield fi;
  uint i, k, n;
  real e, emax;
  
  n = t->n;
  emax = 0.0;
  
  for (i=0; i<n; i++) {
    f = t->mol[i]->f;
    set_zero (3, f);
    k = 3 * i;
    e = 0.0;
    
    fi = flj->v + k;
    axpy (3, 1.0, fi, 1, f, 1);
    e += dot (false, false, 3, fi, 1, fi, 1);
    
    fi = fc->v + 3 * k;
    axpy (3, 1.0, fi, 1, f, 1);
    e += dot (false, false, 3, fi, 1, fi, 1);
    fi = fi + 3;
    axpy (3, 1.0, fi, 1, f, 1);
    e += dot (false, false, 3, fi, 1, fi, 1);
    fi = fi + 3;
    axpy (3, 1.0, fi, 1, f, 1);
    e += dot (false, false, 3, fi, 1, fi, 1);
    
    e = nrm2 (3, f, 1) / REAL_SQRT (e);
    emax = (e > emax) ? e : emax;
  }
  
  return emax;
}

// z <- z + v x w
static void 
cross (real *norm2, pcreal v, pcreal w, preal z) {
  real a, b, c;
  
  a = v[1] * w[2] - v[2] * w[1];
  b = v[2] * w[0] - v[0] * w[2];
  c = v[0] * w[1] - v[1] * w[0];
  
  *norm2 += a * a + b * b + c * c;
  
  z[0] += a;
  z[1] += b;
  z[2] += c;
}

real
calcTorque_tip4p (pavector fc, pavector flj, ptip4p t) {
  pmolecule mol;
  uint i, k, n;
  preal r, rm, com, *xc, tau;
  real e, emax;
  
  n = t->n;
  rm = allocreal (3);
  xc = t->xc;
  emax = 0.0;
  for (i=0; i<n; i++) {
    mol = t->mol[i];
    com = mol->com;
    k = 3 * i;
    e = 0.0;
    tau = mol->t;
    set_zero (3, tau);
    
    r = mol->r[0];
    cross (&e, r, flj->v + k, tau);
    
    // Relative coordinates of the virtual charged site
    rm[0] = xc[k][0] - com[0];
    rm[1] = xc[k][1] - com[1];
    rm[2] = xc[k][2] - com[2];
    cross (&e, rm, fc->v + 3 * k, tau);
    
    r = mol->r[1];
    cross (&e, r, fc->v + 3 * (k + 1), tau);
    
    r = mol->r[2];
    cross (&e, r, fc->v + 3 * (k + 2), tau);
    
    e = nrm2 (3, tau, 1) / REAL_SQRT (e);
    emax = (e > emax) ? e : emax;
  }
  freemem (rm);
  
  return emax;
}

void
update_tip4p (pspatialgeometry sg, pkernelmatrix kc, pkernelmatrix klj, 
              ph2matrix Fc, ph2matrix Flj, pavector fc, pavector flj, ptip4p t) {
  uint i, j, k, l, n;
  preal *xn, *xc, r, com;
  pmolecule *mol = t->mol;
  real d[3], bmin[3], bmax[3], z;
  
  n = t->n;
  xn = t->xn;
  xc = t->xc;
  for (l=0; l<3; l++) {
    bmin[l] = sg->bmin[l];
    bmax[l] = sg->bmax[l];
    d[l] = bmax[l] - bmin[l];
  }
  
  // Recalculate positions of centers of mass, atoms and virtual charges, respectively
  for (i=0; i<n; i++) {
    k = 3 * i;
    com = mol[i]->com;
    r = mol[i]->r[0];
    for (j=0; j<3; j++) {
      xn[i][j] = com[j] + r[j];
    }
    
    r = mol[i]->r[1];
    for (j=0; j<3; j++) {
      xc[k+1][j] = com[j] + r[j];
    }
    
    r = mol[i]->r[2];
    for (j=0; j<3; j++) {
      xc[k+2][j] = com[j] + r[j];
    }
    
    calc_virtual_charge (xc[k+1], xc[k+2], xn[i], xc[k]);
    
    // Check for need of periodic repositioning due to crossing of the unit cell boundary
    for (j=0; j<3; j++) {
      if (com[j] < bmin[j]) {
        z = com[j];
        com[j] += ceil ((bmin[j] - com[j]) / d[j]) * d[j];
        printf ("Moving %le to %le in [%le,%le]\n", z, com[j], bmin[j], bmax[j]);
      }
      else if (com[j] > bmax[j]) {
        z = com[j];
        com[j] -= ceil ((com[j] - bmax[j]) / d[j]) * d[j];
        printf ("Moving %le to %le in [%le,%le]\n", z, com[j], bmin[j], bmax[j]);
      }
      if (xn[i][j] < bmin[j]) {
        z = xn[i][j];
        xn[i][j] += ceil ((bmin[j] - xn[i][j]) / d[j]) * d[j];
        printf ("Moving %le to %le in [%le,%le]\n", z, xn[i][j], bmin[j], bmax[j]);
      }
      else if (xn[i][j] > bmax[j]) {
        z = xn[i][j];
        xn[i][j] -= ceil ((xn[i][j] - bmax[j]) / d[j]) * d[j];
        printf ("Moving %le to %le in [%le,%le]\n", z, xn[i][j], bmin[j], bmax[j]);
      }
    }
    for (l=0; l<3; l++) {
      for (j=0; j<3; j++) {
        if (xc[k+l][j] < bmin[j]) {
          z = xc[k+l][j];
          xc[k+l][j] += ceil ((bmin[j] - xc[k+l][j]) / d[j]) * d[j];
          printf ("Moving %le to %le in [%le,%le]\n", z, xc[k+l][j], bmin[j], bmax[j]);
        }
        else if (xc[k+l][j] > bmax[j]) {
          z = xc[k+l][j];
          xc[k+l][j] -= ceil ((xc[k+l][j] - bmax[j]) / d[j]) * d[j];
          printf ("Moving %le to %le in [%le,%le]\n", z, xc[k+l][j], bmin[j], bmax[j]);
        }
      }
    }
  }
  
  // Recalculate forces and torques
  update_kernelmatrix (3 * n, xc, t->xcmol, kc);
  update_kernelmatrix (n, xn, t->xnmol, klj);
  fullForce_coulomb (t->q, kc, sg, Fc, fc);
  fullForce_lj (klj, sg, Flj, flj);
  calcForce_tip4p (fc, flj, t);
  calcTorque_tip4p (fc, flj, t);
}

/** ---------------------------------------------------------------------
 *    Optimization
 *  --------------------------------------------------------------------- */

real
gradientDescent_tip4p (pspatialgeometry sg, pkernelmatrix kc, pkernelmatrix klj, 
                       ph2matrix Vc, ph2matrix Fc, ph2matrix Vlj, ph2matrix Flj, 
                       real lambda, real nu, uint m, ptip4p t) {
  real E0, E, norm2, b;
  preal *xc_new, *xn_new, *xc, *xn, f;
  uint i, j, k, n, *xcmol, *xnmol;
  pavector q, fc, flj;
  
  n = t->n;
  xn = t->xn;
  xc = t->xc;
  xcmol = t->xcmol;
  xnmol = t->xnmol;
  q = t->q;
  fc = new_avector (9 * n);
  flj = new_avector (3 * n);
  
  // Compute current potential energy and forces
  update_kernelmatrix (3 * n, xc, xcmol, kc);
  update_kernelmatrix (n, xn, xnmol, klj);
  E = E0 = energy_coulomb (q, kc, sg, Vc) + energy_lj (klj, sg, Vlj);
  fullForce_coulomb (q, kc, sg, Fc, fc);
  
  fullForce_lj (klj, sg, Flj, flj);
  calcForce_tip4p (fc, flj, t);
  
  // Set up line search
  norm2 = 0.0;
  for (i=0; i<n; i++) {
    f = t->mol[i]->f;
    norm2 += dot (false, false, 3, f, 1, f, 1);
  }
  b = 0.5 * norm2;
  k = 0;
  xc_new = (preal *) allocmem (3 * n * sizeof (preal));
  xn_new = (preal *) allocmem (n * sizeof (preal));
  xc_new[0] = allocreal (9 * n);
  xn_new[0] = allocreal (3 * n);
  for (i=0; i<n; i++) {
    j = 3 * i;
    xc_new[j] = xc_new[0] + 3 * j;
    xc_new[j+1] = xc_new[j] + 3;
    xc_new[j+2] = xc_new[j] + 6;
    xn_new[i] = xn_new[0] + j;
  }
  
  // Perform line search in direction of the negative gradient
  while (E0 - E < lambda * b && k < m) {
    // Reset coordinates
    for (i=0; i<9*n; i++) {
      xc_new[0][i] = xc[0][i];
    }
    for (i=0; i<3*n; i++) {
      xn_new[0][i] = xn[0][i];
    }
    
    //Compute new coordinates
    for (i=0; i<n; i++) {
      j = 3 * i;
      f = t->mol[i]->f;
      axpy (3, lambda, f, 1, xn_new[i], 1);
      axpy (3, lambda, f, 1, xc_new[j], 1);
      axpy (3, lambda, f, 1, xc_new[j+1], 1);
      axpy (3, lambda, f, 1, xc_new[j+2], 1);
    }
    
    //Compute new potential energy
    update_kernelmatrix (3 * n, xc_new, xcmol, kc);
    update_kernelmatrix (n, xn_new, xnmol, klj);
    E = energy_coulomb (q, kc, sg, Vc) + energy_lj (klj, sg, Vlj);
    lambda *= nu;
    k++;
  }
  
  // Check whether line search was successful and in case update coordinates
  if (E0 - E >= lambda * b) {
    freemem (xc[0]);
    freemem (xn[0]);
    for (i=0; i<n; i++) {
      j = 3 * i;
      xc[j] = xc_new[j];
      xc[j+1] = xc_new[j+1];
      xc[j+2] = xc_new[j+2];
      xn[i] = xn_new[i];
    }
    freemem (xn_new);
    freemem (xc_new);
  }
  else {
    freemem (xc_new[0]);
    freemem (xn_new[0]);
    freemem (xc_new);
    freemem (xn_new);
  }
  del_avector (flj);
  del_avector (fc);

  return E;
}

/*
static void
quaternion_velocity (pcquaternion q, pcreal w, pmolecule mol) {
  pquaternion qv = mol->qv;
  
  qv[0] = -0.5 * (q[1] * w[0] + q[2] * w[1] + q[3] * w[2]);
  qv[1] = 0.5 * (q[0] * w[0] - q[3] * w[1] + q[2] * w[2]);
  qv[2] = 0.5 * (q[3] * w[0] + q[0] * w[1] - q[1] * w[2]);
  qv[3] = 0.5 * (-q[2] * w[0] + q[1] * w[1] + q[0] * w[2]);
}

real
gradientDescentRotation_tip4p (pspatialgeometry sg, pkernelmatrix kc, pkernelmatrix klj, 
                              ph2matrix Vc, ph2matrix Fc, ph2matrix Vlj, ph2matrix Flj, 
                              real lambda, real nu, uint m, ptip4p t) {
  real E0, E, norm2, b, et;
  preal *xc_new, *xn_new, *xc, *xn;
  pfield tau, com, r;
  uint i, j, k, n, *xcmol, *xnmol;
  pavector q, fc, flj;
  pquaternion *p;
  pmolecule *mol = t->mol;
  
  n = t->n;
  xn = t->xn;
  xc = t->xc;
  xcmol = t->xcmol;
  xnmol = t->xnmol;
  q = t->q;
  fc = new_avector (9 * n);
  flj = new_avector (3 * n);
  
  // Compute current potential energy and torques
  update_kernelmatrix (3 * n, xc, xcmol, kc);
  update_kernelmatrix (n, xn, xnmol, klj);
  E = E0 = energy_coulomb (q, kc, sg, Vc) + energy_lj (klj, sg, Vlj);
  fullForce_coulomb (q, kc, sg, Fc, fc);
  
  fullForce_lj (klj, sg, Flj, flj);
  calcTorque_tip4p (fc, flj, t);
  inertia_tip4p (t);
  
  // Set up line search
  norm2 = 0.0;
  
  for (i=0; i<n; i++) {
    refOrientation_molecule (mol[i]);
    tau = mol[i]->t->v;
    norm2 += dot (false, false, 3, tau, 1, tau, 1);
  }
  b = 0.5 * norm2;
  k = 0;
  xc_new = (preal *) allocmem (3 * n * sizeof (preal));
  xn_new = (preal *) allocmem (n * sizeof (preal));
  xc_new[0] = allocreal (9 * n);
  xn_new[0] = allocreal (3 * n);
  p = (pquaternion *) allocmem (n * sizeof (pquaternion));
  for (i=0; i<n; i++) {
    j = 3 * i;
    xc_new[j] = xc_new[0] + 3 * j;
    xc_new[j+1] = xc_new[j] + 3;
    xc_new[j+2] = xc_new[j] + 6;
    xn_new[i] = xn_new[0] + j;
    p[i] = new_quaternion ();
  }
  
  // Perform line search in direction given by the torque
  while (E0 - E < lambda * b && k < m) {
    // Reset quaternions
    for (i=0; i<n; i++) {
      p[i][0] = mol[i]->q[0];
      p[i][1] = mol[i]->q[1];
      p[i][2] = mol[i]->q[2];
      p[i][3] = mol[i]->q[3];
    }
    
    //Compute new coordinates
    for (i=0; i<n; i++) {
      j = 3 * i;
      tau = mol[i]->t->v;
      com = mol[i]->com->v;
      quaternion_velocity (p[i], tau, mol[i]);
      axpy (4, lambda, mol[i]->qv, 1, p[i], 1);
      normalize_quaternion (p[i]);
      adjust_molecule (p[i], mol[i]);
      r = mol[i]->r[0]->v;
      xn_new[i][0] = com[0] + r[0];
      xn_new[i][1] = com[1] + r[1];
      xn_new[i][2] = com[2] + r[2];
      r = mol[i]->r[1]->v;
      xc_new[j+1][0] = com[0] + r[0];
      xc_new[j+1][1] = com[1] + r[1];
      xc_new[j+1][2] = com[2] + r[2];
      r = mol[i]->r[2]->v;
      xc_new[j+2][0] = com[0] + r[0];
      xc_new[j+2][1] = com[1] + r[1];
      xc_new[j+2][2] = com[2] + r[2];
      calc_virtual_charge (xc_new[j+1], xc_new[j+2], xn_new[i], xc_new[j]);
    }
    
    //Compute new potential energy
    update_kernelmatrix (3 * n, xc_new, xcmol, kc);
    update_kernelmatrix (n, xn_new, xnmol, klj);
    E = energy_coulomb (q, kc, sg, Vc) + energy_lj (klj, sg, Vlj);
    
    lambda *= nu;
    k++;
  }
  
  // Check whether line search was successful and in case update coordinates
  if (E0 - E >= lambda * b) {
    freemem (xc[0]);
    freemem (xn[0]);
    // Accept new coordinates
    for (i=0; i<n; i++) {
      j = 3 * i;
      xc[j] = xc_new[j];
      xc[j+1] = xc_new[j+1];
      xc[j+2] = xc_new[j+2];
      xn[i] = xn_new[i];
      mol[i]->q[0] = p[i][0];
      mol[i]->q[1] = p[i][1];
      mol[i]->q[2] = p[i][2];
      mol[i]->q[3] = p[i][3];
      del_quaternion (p[i]);
    }
    freemem (xn_new);
    freemem (xc_new);
    freemem (p);
  }
  else {
    // Reject new coordinates
    for (i=0; i<n; i++) {
      del_quaternion (p[i]);
      adjust_molecule (mol[i]->q, mol[i]);
    }
    freemem (xc_new[0]);
    freemem (xn_new[0]);
    freemem (xc_new);
    freemem (xn_new);
    freemem (p);
  }
  
  del_avector (flj);
  del_avector (fc);

  return E;
}

*/

/** ----------------------------------------------------------------------
 *    Time-stepping methods
 *  ---------------------------------------------------------------------- */

void
initMolecules_tip4p (ptip4p t) {
  uint i, n;
  pmolecule mol;
  preal I, *r;
  real r00, r01, r02, r10, r11, r12, r20, r21, r22;
  
  n = t->n;
  for (i=0; i<n; i++) {
    mol = t->mol[i];
    I = mol->R;
    r = mol->r;
    r00 = r[0][0];
    r01 = r[0][1];
    r02 = r[0][2];
    r10 = r[1][0];
    r11 = r[1][1];
    r12 = r[1][2];
    r20 = r[2][0];
    r21 = r[2][1];
    r22 = r[2][2];
    
    // Diagonal entries
    I[0] = M_O * (r01 * r01 + r02 * r02) 
                    + M_H * (r11 * r11 + r12 * r12 + r21 * r21 + r22 * r22);
    I[4] = M_O * (r00 * r00 + r02 * r02) 
                    + M_H * (r10 * r10 + r12 * r12 + r20 * r20 + r22 * r22);
    I[8] = M_O * (r00 * r00 + r01 * r01) 
                    + M_H * (r10 * r10 + r11 * r11 + r20 * r20 + r21 * r21);
                  
    // Off-diagonal entries
    I[1] = I[3] = - M_O * r00 * r01 - M_H * (r10 * r11 + r20 * r21);
    I[2] = I[6] = - M_O * r00 * r02 - M_H * (r10 * r12 + r20 * r22);
    I[5] = I[7] = - M_O * r01 * r02 - M_H * (r11 * r12 + r21 * r22);
    
    // Reference orientation
    refOrientation_molecule (mol);
  }
}

real
velocityVerlet_tip4p (pspatialgeometry sg, pkernelmatrix kc, pkernelmatrix klj, 
                      ph2matrix Fc, ph2matrix Flj, pavector fc, pavector flj, 
                      real delta, ptip4p t) {
  uint n = t->n;
  pmolecule *mol = t->mol;
  uint i;
  real E;
  
  for (i=0; i<n; i++) {
    positionVerlet_molecule (delta, mol[i]);
  }
  update_tip4p (sg, kc, klj, Fc, Flj, fc, flj, t);
  E = 0.0;
  for (i=0; i<n; i++) {
    E += velocityVerlet_molecule (delta, mol[i]);
  }
  
  return E;
}