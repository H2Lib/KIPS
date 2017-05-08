
/* ------------------------------------------------------------
 * This is the file "particles.c" of the KIPS package.
 * All rights reserved, Steffen Boerm 2017
 * ------------------------------------------------------------ */

#include "particles.h"
#include "basic.h"
#include "interpolation.h"

pparticles
new_particles(uint n, uint m)
{
  pparticles p;
  real **x;
  real *xdata;
  uint i;

  p = (pparticles) allocmem(sizeof(particles));
  p->n = n;
  p->x = x = (real **) allocmem(sizeof(real *) * n);
  p->xdata = xdata = allocreal(n * 3);

  for(i=0; i<n; i++)
    p->x[i] = xdata + 3*i;

  p->in = new_interpolation(m);
  chebyshev_interpolation(-1.0, 1.0, p->in);
  
  return p;
}

void
del_particles(pparticles p)
{
  del_interpolation(p->in);
  p->in = 0;

  freemem(p->xdata);
  p->xdata = 0;

  freemem(p->x);
  p->x = 0;

  freemem(p);
}

void
random_particles(pparticles p)
{
  real **x = p->x;
  uint n = p->n;
  uint i;

  for(i=0; i<n; i++) {
    x[i][0] = REAL_RAND();
    x[i][1] = REAL_RAND();
    x[i][2] = REAL_RAND();
  }
}

pclustergeometry
buildgeometry_particles(pcparticles p)
{
  pclustergeometry cg;
  const real **x = (const real **) p->x;
  uint n = p->n;
  uint i;

  cg = new_clustergeometry(3, n);

  for(i=0; i<n; i++) {
    cg->x[i][0] = cg->smin[i][0] = cg->smax[i][0] = x[i][0];
    cg->x[i][1] = cg->smin[i][1] = cg->smax[i][1] = x[i][1];
    cg->x[i][2] = cg->smin[i][2] = cg->smax[i][2] = x[i][2];
  }

  return cg;
}

field
potential_newton(const real *x, const real *y)
{
  return REAL_RSQRT(REAL_SQR(x[0]-y[0]) +
		    REAL_SQR(x[1]-y[1]) +
		    REAL_SQR(x[2]-y[2]));
}

void
addeval_direct_particles(field alpha, pcparticles p,
			 pcavector x, pavector y)
{
  const real **px = (const real **) p->x;
  uint n = p->n;
  pcfield vx = x->v;
  pfield vy = y->v;
  uint i, j;

  for(j=0; j<n; j++)
    for(i=0; i<n; i++)
      vy[i] += alpha * potential_newton(px[i], px[j]) * vx[j];
}

void
buildV_particles(pcparticles p, pccluster t, pamatrix V)
{
  const real **x = (const real **) p->x;
  pcinterpolation in = p->in;
  uint m = p->in->m;
  pfield Va;
  longindex ldV;
  pinterpolation in0, in1, in2;
  uint nu0, nu1, nu2;
  real lagr, lagr1, lagr2;
  uint i, ii, j, j1, j2;

  /* Adjust size of V if necessary */
  if(V->rows != t->size || V->cols != m*m*m)
    resize_amatrix(V, t->size, m*m*m);

  Va = V->a;
  ldV = V->ld;

  /* Set up interpolation points for x, y, and z coordinates */
  in0 = build_transformed_interpolation(in, t->bmin[0], t->bmax[0]);
  in1 = build_transformed_interpolation(in, t->bmin[1], t->bmax[1]);
  in2 = build_transformed_interpolation(in, t->bmin[2], t->bmax[2]);

  /* Iterate over all cluster indices */
  for(i=0; i<t->size; i++) {
    /* Get particle index */
    ii = t->idx[i];

    for(nu2=0; nu2<m; nu2++) {
      /* Evaluate one-dimensional Lagrange polynomial in z */
      lagr2 = evallagrange_interpolation(in2, nu2, x[ii][2]);
      j2 = nu2;

      for(nu1=0; nu1<m; nu1++) {
	/* Evaluate one-dimensional Lagrange polynomial in y */
	lagr1 = lagr2 * evallagrange_interpolation(in1, nu1, x[ii][1]);
	j1 = nu1 + m * j2;

	for(nu0=0; nu0<m; nu0++) {
	  /* Evaluate one-dimensional Lagrange polynomial in x */
	  lagr = lagr1 * evallagrange_interpolation(in0, nu0, x[ii][0]);
	  j = nu0 + m * j1;

	  Va[i+j*ldV] = lagr;
	}
      }
    }
  }

  /* Clean up interpolation points */
  del_interpolation(in2);
  del_interpolation(in1);
  del_interpolation(in0);
}

void
buildE_particles(pcparticles p, pccluster son, pccluster father,
		 pamatrix E)
{
  pcinterpolation in = p->in;
  uint m = p->in->m;
  pfield Ea;
  longindex ldE;
  pinterpolation sin0, sin1, sin2;
  pinterpolation fin0, fin1, fin2;
  uint nu0, nu1, nu2;
  uint mu0, mu1, mu2;
  real lagr, lagr1, lagr2;
  uint i, i1, i2, j, j1, j2;

  /* Adjust size of V if necessary */
  if(E->rows != m*m*m || E->cols != m*m*m)
    resize_amatrix(E, m*m*m, m*m*m);

  Ea = E->a;
  ldE = E->ld;

  /* Set up interpolation points for x, y, and z coordinates
   * of the son cluster */
  sin0 = build_transformed_interpolation(in, son->bmin[0], son->bmax[0]);
  sin1 = build_transformed_interpolation(in, son->bmin[1], son->bmax[1]);
  sin2 = build_transformed_interpolation(in, son->bmin[2], son->bmax[2]);

  /* Set up interpolation points for x, y, and z coordinates
   * of the son cluster */
  fin0 = build_transformed_interpolation(in, father->bmin[0], father->bmax[0]);
  fin1 = build_transformed_interpolation(in, father->bmin[1], father->bmax[1]);
  fin2 = build_transformed_interpolation(in, father->bmin[2], father->bmax[2]);

  /* Iterate over z coordinates */
  for(mu2=0; mu2<m; mu2++)
    for(nu2=0; nu2<m; nu2++) {
      j2 = mu2;
      i2 = nu2;
      lagr2 = evallagrange_interpolation(fin2, mu2, sin2->xi[nu2]);

      /* Iterate over y coordinates */
      for(mu1=0; mu1<m; mu1++)
	for(nu1=0; nu1<m; nu1++) {
	  j1 = mu1 + m * j2;
	  i1 = nu1 + m * i2;
	  lagr1 = lagr2 * evallagrange_interpolation(fin1, mu1, sin1->xi[nu1]);

	  /* Iterate over x coordinates */
	  for(mu0=0; mu0<m; mu0++)
	    for(nu0=0; nu0<m; nu0++) {
	      j = mu0 + m * j1;
	      i = nu0 + m * i1;
	      lagr = lagr1 * evallagrange_interpolation(fin0, mu0, sin0->xi[nu0]);

	      Ea[i+j*ldE] = lagr;
	    }
	}
    }

  /* Clean up interpolation points */
  del_interpolation(fin2);
  del_interpolation(fin1);
  del_interpolation(fin0);
  del_interpolation(sin2);
  del_interpolation(sin1);
  del_interpolation(sin0);
}
