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
  freemem (sg->bmin);
  freemem (sg->bmax);
  freemem (sg->splits);
  freemem (sg->s);
  freemem (sg);
}

/* ------------------------------------------------------------
   Spatial cluster tree
   ------------------------------------------------------------ */

pspatialcluster
init_spatialgeometry(uint maxdepth, real maxdiam, pspatialgeometry sg)
{
  preal bmin, bmax, _bmin, d;
  uint i, k, depth, l, m, r, s, dim, desc;
  uint *n, *j, *n_prod;
  real diam, x, y;
  uint** splits_temp;
  bool fw;
  
  assert (maxdepth > 0 && maxdiam > 0.0);
  
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
  if (diam > maxdiam) {
    depth = dim * ceil (REAL_LOG(diam/maxdiam)/REAL_LOG(2.0));
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
  while (depth < maxdepth && diam > maxdiam) {
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
  for (l=0;l<=depth;l++) {
    sg->splits[l] = allocuint(dim);
    for (i=0; i<dim; i++) {
      sg->splits[l][i] = splits_temp[l][i];
    }
  }
  
  /* Set up spatial clusters */
  d = allocreal (dim);
  n = allocuint (dim);
  n_prod = allocuint (dim+1);
  j = allocuint (dim);
  _bmin = allocreal (dim);
  for (i=0; i<dim; i++) {
    _bmin[i] = bmin[i] = sg->bmin [i];
    bmax[i] = sg->bmax[i];
    d[i] = bmax[i]-bmin[i];
    n[i] = 1;
    j[i] = 0;
    n_prod[i] = 1;
  }
  n_prod[dim] = 1;
  desc = (1<<(depth+1)) - 1;
  sg->s = (pspatialcluster **) allocmem((size_t) (depth+1)*sizeof(pspatialcluster *));
  for (l=0; l<=depth; l++) {
    if (l>=1) {
      for (i=0; i<dim; i++) {
        if (splits_temp[l][i]-splits_temp[l-1][i]==1) {
          m = i;
        }
      }
      n[m] *= 2;
      d[m] /= 2.0;
      bmax[m] = bmin[m]+d[m];
      for (i=m+1; i<=dim; i++) {
        n_prod[i] *= 2;
      }
      desc -= 1<<(depth-l+1);
    }
    
    sg->s[l] = (pspatialcluster *) allocmem((size_t) n_prod[dim]*sizeof(pspatialcluster));
    
    k = 0;
    if (l == depth) {
        sg->s[l][k] = new_spatialcluster (dim, bmin, bmax,0,desc);
    }
    else {
      sg->s[l][k] = new_spatialcluster (dim, bmin, bmax,2,desc);
    }
    
    if (l>=1) {
          sg->s[l-1][k]->son[0] = sg->s[l][k];
    }
    
    k++;
    fw = true;;
    while (fw) {
      for (i=0; i<dim; i++) {
        j[i]++;
        if (j[i] < n[i]) {
          bmin[i] += d[i];
          bmax[i] += d[i];
          if (l == depth) {
            sg->s[l][k] = new_spatialcluster (dim, bmin, bmax,0,desc);
          }
          else {
            sg->s[l][k] = new_spatialcluster (dim, bmin, bmax,2,desc);
          }
          
          if (l>=1) {
            r = 0;
              for (s=m+1; s<dim; s++) {
                r += j[s] * n_prod[s];
              }
              r*= 0.5;
              if (j[m]%2 == 0) {
                sg->s[l-1][k-r-j[m]/2*n_prod[m]]->son[0] = sg->s[l][k];
              }
              else {
                sg->s[l-1][k-r-(j[m]+1)/2*n_prod[m]]->son[1] = sg->s[l][k];
              }
              
          }
            
          k++;
          fw = true;
          break;
        }
        else {
          j[i] = 0;
          bmin[i] = _bmin[i];
          bmax[i] = bmin[i]+d[i];
          fw = false;
        }
      }
    }
  }
  freemem (n);
  freemem (j);
  freemem (bmin);
  freemem (bmax);
  freemem (_bmin);
  freemem (n_prod);
  freemem (d);
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
  
  if (nidx >= 1) {
    depth = sg->depth;
    n = countonlevel_spatialgeometry (depth, sg);
    
    idx = allocuint (nidx);
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
