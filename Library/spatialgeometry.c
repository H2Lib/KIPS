/* ------------------------------------------------------------
 * This is the file "spatialgeometry.c" of the KIPS package.
 * All rights reserved, Jonas Lorenzen 2018
 * ------------------------------------------------------------ */

#include "spatialgeometry.h"
#include "basic.h"
#include <assert.h>
#include <math.h>


/* ------------------------------------------------------------
 * Constructors and destructors
 * ------------------------------------------------------------ */

pspatialgeometry
new_spatialgeometry(uint dim, preal bmin, preal bmax)
{
  pspatialgeometry sg;
  
  sg = (pspatialgeometry) allocmem(sizeof(spatialgeometry));
  sg->bmin = bmin;
  sg->bmax = bmax;
  sg->dim = dim;
  
  return sg;
}

void
del_spatialgeometry(pspatialgeometry sg)
{
  uint l;
  freemem (sg->bmin);
  freemem (sg->bmax);
  freemem (sg->s[0][0]->idx);
  for (l=0; l<=sg->depth; l++) {
    freemem (sg->splits[l]);
    freemem (sg->s[l]);
  }
  freemem (sg->splits);
  freemem (sg->s);
  freemem (sg);
}

/* ------------------------------------------------------------
   Spatial cluster tree
   ------------------------------------------------------------ */

/* Enumerate clusters on a given level recursively by dimension and 
   set up bounding boxes and links to sons. */
static void
label_cluster (uint level, uint dim, uint m, uint k, uint k_son0, 
               uint k_son1, pspatialgeometry sg) {
  uint i, n, sons;
  real d;
  pspatialcluster s;
  
  assert (level >= 1);
  dim --;
  n = 1<<sg->splits[level-1][dim];
  sons = (level == sg->depth ? 0 : 2);
  if (dim < sg->dim) {
    k *= n;
    k_son0 *= n;
    k_son1 *= n;
  }
  if (dim > 0) {
    if (dim == m) {
      for (i=0; i<n; i++) {
        label_cluster (level, dim, m, i+k, 2*(i+k_son0), 2*(i+k_son1)+1, sg);
      }
    }
    else {
      for (i=0; i<n; i++) {
        label_cluster (level, dim, m, i+k, i+k_son0, i+k_son1, sg);
      }
    }
  }
  else {
    if (dim == m) {
      for (i=0; i<n; i++) {
        s = sg->s[level-1][i+k];
        s->son[0] = sg->s[level][2*(i+k_son0)] 
                  = new_spatialcluster (sg->dim, s->bmin, s->bmax, sons);
        s->son[1] = sg->s[level][2*(i+k_son1)+1]
                  = new_spatialcluster (sg->dim, s->bmin, s->bmax, sons);
        d = 0.5 * (s->bmax[m]-s->bmin[m]);
        s->son[0]->bmax[m] -= d;
        s->son[1]->bmin[m] += d; 
      }
    }
    else {
      for (i=0; i<n; i++) {
        s = sg->s[level-1][i+k];
        s->son[0] = sg->s[level][i+k_son0] 
                  = new_spatialcluster (sg->dim, s->bmin, s->bmax, sons);
        s->son[1] = sg->s[level][i+k_son1]
                  = new_spatialcluster (sg->dim, s->bmin, s->bmax, sons);
        d = 0.5 * (s->bmax[m]-s->bmin[m]);
        s->son[0]->bmax[m] -= d;
        s->son[1]->bmin[m] += d; 
      }
    }
  }
}

pspatialcluster
init_spatialgeometry(uint maxdepth, real mindiam, pspatialgeometry sg)
{
  preal bmin, bmax;
  uint i, k, depth, level, m, dim;
  real diam, x, y;
  uint** splits_temp;
  
  assert (maxdepth > 0 && mindiam > 0.0);
  
  dim = sg->dim;
  diam = 0.0;
  bmin = allocreal(dim);
  bmax = allocreal(dim);
  for (i=0; i<dim; i++) {
    bmin[i] = sg->bmin[i];
    bmax[i] = sg->bmax[i];
    diam += REAL_SQR(bmax[i]-bmin[i]);
  }
  diam = REAL_SQRT(diam);
  
  /* Find upper bound for the needed depth. */
  if (diam > mindiam) {
    depth = dim * ceil (REAL_LOG(diam/mindiam)/REAL_LOG(2.0));
  }
  else {
    depth = 0;
  }
  if (depth > maxdepth) {
    depth = maxdepth;
  }
  splits_temp = allocmem((size_t) (depth+1)*sizeof(uint *));
  
  /* Determine the actual depth and needed bisection steps. */
  depth = 0;
  splits_temp[0] = allocuint(dim);
  for (i=0; i<dim; i++) {
    splits_temp[0][i] = 0;
  }
  while (depth < maxdepth && diam > mindiam) {
    x = 0.0;
    for (i=0; i<dim; i++) {
      y = bmax[i]-bmin[i];
      if (y > x) {
        x = y;
        k = i;
      }
    }
    
    bmax[k] -= 0.5*x;
    depth ++;
    splits_temp[depth] = allocuint(dim);
    for (i=0; i<dim; i++) {
      if (i==k) {
        splits_temp[depth][i] = splits_temp[depth-1][i]+1;
      }
      else {
        splits_temp[depth][i] = splits_temp[depth-1][i];
      }
    }
    diam = REAL_SQRT(REAL_SQR(diam)-0.75*REAL_SQR(x));
  }
  sg->depth = depth;
  sg->splits = (uint **) allocmem((size_t) (depth+1)*sizeof(uint *));
  for (level=0; level<=depth; level++) {
    sg->splits[level] = allocuint(dim);
    for (i=0; i<dim; i++) {
      sg->splits[level][i] = splits_temp[level][i];
    }
  }
  
  /* Set up spatial clusters */
  
  sg->s = (pspatialcluster **) allocmem((size_t) (depth+1)*sizeof(pspatialcluster *));
  sg->s[0] = (pspatialcluster *) allocmem ((size_t) (sizeof (pspatialcluster)));
  sg->s[0][0] = (depth == 0 ? new_spatialcluster (dim, sg->bmin, sg->bmax, 0) : 
                 new_spatialcluster (dim, sg->bmin, sg->bmax, 2));
  for (level=1; level<=depth; level++) {
    sg->s[level] 
      = (pspatialcluster *) allocmem ((size_t) ((1<<level) * sizeof (pspatialcluster)));
    for (i=0; i<dim; i++) {
      if (splits_temp[level][i]-splits_temp[level-1][i]==1) {
        m = i;
      }
    }
    label_cluster (level, dim, m, 0, 0, 0, sg);
  }
  
  update_spatialcluster (sg->s[0][0]);
  
  freemem (bmin);
  freemem (bmax);
  for (level=0; level<=depth; level++) {
    freemem (splits_temp[level]);
  }
  freemem (splits_temp);
  
  return sg->s[0][0];
}

pspatialcluster
findCluster_spatialgeometry(uint l, pcreal x, pcspatialgeometry sg, uint *nr)
{
  uint i, j, k, n, n_prod, dim;
  real min, max, d;
  
  assert (l <= sg->depth);
  
  dim = sg->dim;
  k = 0;
  n = n_prod = 1;
  for (i=0; i<dim; i++) {
    n_prod *= n;
    min = sg->bmin[i];
    max = sg->bmax[i];
    d = max-min;
    assert (d > 0.0);
    assert (x[i] >= min && x[i] < max);
    n = countindirection_spatialgeometry (l, i, sg);
    
    j = floor ((x[i]-min)/d*n);
    k += j*n_prod;
  }
  
  if (nr) {
    *nr = k;
  }
  return sg->s[l][k];
}

void
initPoints_spatialgeometry (uint nidx, pcreal *x, pspatialgeometry sg) {
  uint *idx, *count;
  uint i, j, n, nr, depth;
  pspatialcluster s;
  
  freemem (sg->s[0][0]->idx);
  if (nidx >= 1) {
    depth = sg->depth;
    n = countonlevel_spatialgeometry (depth, sg);
    
    sg->s[0][0]->idx = allocuint (nidx);
    idx = sg->s[0][0]->idx;
    count = allocuint (n);
    for (j=0; j<n; j++) {
      count[j] = 0;
    }
    for (i=0; i<nidx; i++) {
      s = findCluster_spatialgeometry (depth, x[i], sg, &nr);
      count[nr]++;
    }
    sg->s[depth][0]->idx = idx;
    sg->s[depth][0]->nidx = count[0];
    for (j=1; j<n; j++) {
      sg->s[depth][j]->idx = idx + count[j-1];
      sg->s[depth][j]->nidx = count[j];
      count[j] += count[j-1];
    }
    for (i=0; i<nidx; i++) {
      s = findCluster_spatialgeometry (depth, x[i], sg, &nr);
      assert (count[nr] > 0);
      idx[count[nr]-1] = i;
      count[nr]--;
    }
  }
  update_spatialcluster (sg->s[0][0]);
  
  (void) s;
  freemem (count);
}

/* ------------------------------------------------------------
   Derived quantities
   ------------------------------------------------------------ */

uint
countonlevel_spatialgeometry (uint l, pcspatialgeometry sg) {
  uint i, n, dim;
  
  assert (l <= sg->depth);
  
  dim = sg->dim;
  n = 1;
  for (i=0; i<dim; i++) {
    n *= 1<<sg->splits[l][i];
  }
  
  return n;
}

uint
countindirection_spatialgeometry (uint l, uint d, pcspatialgeometry sg) {
  assert (l <= sg->depth);
  assert (d < sg->dim);
  
  return 1<<sg->splits[l][d];
}
