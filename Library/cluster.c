
/* ------------------------------------------------------------
 * This is the file "cluster.c" of the KIPS package.
 * All rights reserved, Steffen Boerm 2010
 * ------------------------------------------------------------ */

#include "cluster.h"

#include "basic.h"

#include <assert.h>
#include <math.h>

#ifdef USE_NETCDF
#include <stdio.h>
#include <string.h>
#include <netcdf.h>
#endif

static uint active_cluster = 0;

/* ------------------------------------------------------------
 * Constructors and destructors
 * ------------------------------------------------------------ */

pcluster
new_cluster(uint size, uint *idx, uint sons, uint dim)
{
  pcluster t;
  uint i;

  t = (pcluster) allocmem(sizeof(cluster));
  t->size = size;
  t->idx = idx;
  t->sons = sons;
  t->dim = dim;
  t->desc = 0;

  t->son = NULL;
  if(sons > 0) {
    t->son = (pcluster *) allocmem((size_t) sizeof(cluster) * sons);
    for(i=0; i<sons; i++)
      t->son[i] = NULL;
  }

  t->bmin = allocreal(dim);
  t->bmax = allocreal(dim);

#ifdef USE_OPENMP
#pragma omp atomic
#endif
  active_cluster++;

  return t;
}

void
update_cluster(pcluster t)
{
  uint desc;
  uint i;

  desc = 1;
  for(i=0; i<t->sons; i++)
    desc += t->son[i]->desc;

  t->desc = desc;
}

void
del_cluster(pcluster t)
{
  uint i;
      
  freemem(t->bmin);
  freemem(t->bmax);

  if(t->sons > 0) {
    for(i=0; i<t->sons; i++)
      del_cluster(t->son[i]);
    freemem(t->son);
  }

  freemem(t);

  assert(active_cluster > 0);

#ifdef USE_OPENMP
#pragma omp atomic
#endif
  active_cluster--;
}

void
setsons_cluster(pcluster t, uint sons)
{
  uint i;

  assert(t->sons == 0);

  t->sons = sons;
  t->son = (pcluster *) allocmem(sizeof(pcluster) * sons);

  for(i=0; i<sons; i++)
    t->son[i] = 0;
}

/* ------------------------------------------------------------
 * Bounding box
 * ------------------------------------------------------------ */

real
diam_cluster(pccluster t)
{
  real diam2;
  uint i;

  diam2 = 0.0;
  for(i=0; i<t->dim; i++)
    diam2 += REAL_SQR(t->bmax[i] - t->bmin[i]);

  return REAL_SQRT(diam2);
}

real
dist_cluster(pccluster t, pccluster s)
{
  real dist2;
  uint i;

  assert(t->dim == s->dim);

  dist2 = 0.0;
  for(i=0; i<t->dim; i++)
    if(t->bmax[i] < s->bmin[i])
      dist2 += REAL_SQR(s->bmin[i] - t->bmax[i]);
    else if(s->bmax[i] < t->bmin[i])
      dist2 += REAL_SQR(t->bmin[i] - s->bmax[i]);

  return REAL_SQRT(dist2);
}

real
getmindiam_cluster(pccluster t)
{
  real mindiam;
  uint i;
  
  if(t->sons > 0) {
    mindiam = getmindiam_cluster(t->son[0]);
    for(i=1; i<t->sons; i++)
      mindiam = REAL_MIN(mindiam, getmindiam_cluster(t->son[i]));
  }
  else {
    mindiam = t->bmax[0] - t->bmin[0];
    for(i=1; i<t->dim; i++)
      mindiam = REAL_MIN(mindiam, t->bmax[i] - t->bmin[i]);
  }

  return mindiam;
}

/* ------------------------------------------------------------
 * Statistics
 * ------------------------------------------------------------ */

uint
getactives_cluster()
{
  return active_cluster;
}

size_t
getsize_cluster(pccluster t)
{
  size_t sz;
  uint i;

  sz = (size_t) sizeof(cluster) + 2 * t->dim * sizeof(real);
  if(t->sons > 0) {
    sz += (size_t) sizeof(pcluster) * t->sons;
    for(i=0; i<t->sons; i++)
      sz += getsize_cluster(t->son[i]);
  }

  return sz;
}

uint
getdepth_cluster(pccluster t)
{
  uint d, d1;
  uint i;

  d = 0;

  if(t->sons) {
    for(i=0; i<t->sons; i++) {
      d1 = getdepth_cluster(t->son[i]);
      if(d1 > d)
	d = d1;
    }
    d++;
  }
  
  return d;
}
