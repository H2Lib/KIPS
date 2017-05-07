
/* ------------------------------------------------------------
 * This is the file "particles.c" of the KIPS package.
 * All rights reserved, Steffen Boerm 2017
 * ------------------------------------------------------------ */

#include "particles.h"
#include "basic.h"

pparticles
new_particles(uint n)
{
  pparticles p;

  p = (pparticles) allocmem(sizeof(particles));
  p->n = n;
  p->x0 = allocreal(n);
  p->x1 = allocreal(n);
  p->x2 = allocreal(n);

  return p;
}

void
del_particles(pparticles p)
{
  freemem(p->x2);
  freemem(p->x1);
  freemem(p->x0);
  freemem(p);
}

void
random_particles(pparticles p)
{
  real *x0 = p->x0;
  real *x1 = p->x1;
  real *x2 = p->x2;
  uint n = p->n;
  uint i;

  for(i=0; i<n; i++) {
    x0[i] = REAL_RAND();
    x1[i] = REAL_RAND();
    x2[i] = REAL_RAND();
  }
}

pclustergeometry
buildgeometry_particles(pcparticles p)
{
  pclustergeometry cg;
  const real *x0 = p->x0;
  const real *x1 = p->x1;
  const real *x2 = p->x2;
  uint n = p->n;
  uint i;

  cg = new_clustergeometry(3, n);

  for(i=0; i<n; i++) {
    cg->x[i][0] = cg->smin[i][0] = cg->smax[i][0] = x0[i];
    cg->x[i][1] = cg->smin[i][1] = cg->smax[i][1] = x1[i];
    cg->x[i][2] = cg->smin[i][2] = cg->smax[i][2] = x2[i];
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
			 pcavector m, pavector phi)
{
  const real *x0 = p->x0;
  const real *x1 = p->x1;
  const real *x2 = p->x2;
  uint n = p->n;
  pcfield vm = m->v;
  pfield vphi = phi->v;
  real x[3], y[3];
  uint i, j;

  for(j=0; j<n; j++) {
    y[0] = x0[j];
    y[1] = x1[j];
    y[2] = x2[j];

    for(i=0; i<n; i++) {
      x[0] = x0[i];
      x[1] = x1[i];
      x[2] = x2[i];

      vphi[i] += alpha * potential_newton(x, y) * vm[j];
    }
  }
}
