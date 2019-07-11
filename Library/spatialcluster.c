/* ------------------------------------------------------------
 * This is the file "spatialcluster.c" of the KIPS package.
 * All rights reserved, Jonas Lorenzen 2018
 * ------------------------------------------------------------ */

#include <assert.h>
#include <stdio.h>

#include "basic.h"
#include "spatialcluster.h"

static uint active_spatialcluster = 0;

/* ------------------------------------------------------------
 * Constructors and destructors
 * ------------------------------------------------------------ */

pspatialcluster
new_spatialcluster(uint dim, preal bmin, preal bmax, uint sons, uint desc) 
{
  pspatialcluster s;
  uint i;
  
  s = (pspatialcluster) allocmem(sizeof(spatialcluster));
  s->dim = dim;
  s->bmin = allocreal(dim);
  s->bmax = allocreal(dim);
  for (i=0; i<dim; i++) {
    s->bmin[i] = bmin[i];
    s->bmax[i] = bmax[i];
  }
  s->idx = NULL;
  s->nidx = 0;
  
  s->sons = sons;
  s->son = NULL;
  if(sons > 0) {
    s->son = (pspatialcluster *) allocmem((size_t) sizeof(spatialcluster) * sons);
    for(i=0; i<sons; i++)
      s->son[i] = NULL;
  }
  s->desc = desc;
  
  active_spatialcluster++;

  return s;
}

void
del_spatialcluster(pspatialcluster s)
{
  uint i;
      
  freemem(s->bmin);
  freemem(s->bmax);
  //freemem(s->idx);

  if(s->sons > 0) {
    for(i=0; i<s->sons; i++)
      del_spatialcluster(s->son[i]);
    freemem(s->son);
  }

  freemem(s);

  assert(active_spatialcluster > 0);

  active_spatialcluster--;
}

void
update_spatialcluster (pspatialcluster s) {
  uint i, nidx;
  
  if (s->sons != 0) {
    nidx = 0;
    for (i=0; i<s->sons; i++) {
      update_spatialcluster (s->son[i]);
      nidx += s->son[i]->nidx;
    }
    s->nidx = nidx;
  }
}

/* ------------------------------------------------------------
   Bounding box
   ------------------------------------------------------------ */

real
diam_spatialcluster(pcspatialcluster s)
{
  real diam2;
  uint i, dim;
  
  dim = s->dim;

  diam2 = 0.0;
  for(i=0; i<dim; i++)
    diam2 += REAL_SQR(s->bmax[i] - s->bmin[i]);

  return REAL_SQRT(diam2);
}

real
dist_spatialcluster(pcspatialcluster s, pcspatialcluster t)
{
  real dist2;
  uint i, dim;

  dim = s->dim;
  assert (dim == t->dim);
  
  dist2 = 0.0;
  for(i=0; i<dim; i++)
    if(t->bmax[i] < s->bmin[i])
      dist2 += REAL_SQR(s->bmin[i] - t->bmax[i]);
    else if(s->bmax[i] < t->bmin[i])
      dist2 += REAL_SQR(t->bmin[i] - s->bmax[i]);

  return REAL_SQRT(dist2);
}

real
distPeriodic_spatialcluster (pcspatialcluster s, pcspatialcluster t, pcreal a, pcreal b) {
  real dist2, r, z;
  uint i, dim;
  preal x, y, v, w;
  
  x = t->bmin;
  y = t->bmax;
  v = s->bmin;
  w = s->bmax;
  
  dim = s->dim;
  assert (dim == t->dim);
  
  dist2 = 0.0;

  for (i=0; i<dim; i++) {
    if (w[i] < x[i]) {
      r = x[i] - w[i];
      z = v[i] - y[i] + b[i] - a[i];
    }
    else if (y[i] < v[i]) {
      r = v[i] - y[i];
      z = x[i] - w[i] + b[i] - a[i];
    }
    else {
      r = z = 0.0;
    }
    dist2 += (r > z) ? REAL_SQR (z) : REAL_SQR (r);
  }
  return REAL_SQRT(dist2);
}

real
distPeriodicPoint_spatialcluster (pcreal x, pcspatialcluster s, pcreal a, pcreal b) {
  real dist2, r, z;
  uint i, dim;
  preal v, w;
  
  v = s->bmin;
  w = s->bmax;
  dim = s->dim;
  
  dist2 = 0.0;
  
  for (i=0; i<dim; i++) {
    assert (a[i] <= x[i] && x[i] <= b[i]);
    if (x[i] < v[i]) {
      r = v[i] - x[i];
      z = x[i] - w[i] + b[i] - a[i];
    }
    else if (x[i] > w[i]) {
      r = x[i] - w[i];
      z = v[i] - x[i] + b[i] - a[i];
    }
    else {
      r = z = 0.0;
    }
    dist2 += (r > z) ? REAL_SQR (z) : REAL_SQR (r);
  }
  return REAL_SQRT(dist2);
}

/* ------------------------------------------------------------
   Statistics
   ------------------------------------------------------------ */

uint
getactives_spatialcluster()
{
  return active_spatialcluster;
}
