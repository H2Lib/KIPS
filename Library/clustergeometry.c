
/* ------------------------------------------------------------
 * This is the file "clustergeometry.h" of the H2Lib package.
 * All rights reserved, Steffen Boerm 2010
 * ------------------------------------------------------------ */

#include "clustergeometry.h"

#include "basic.h"

#include <assert.h>

/* Maximal aspect ratio for bounding boxes */
#define BBOX_MAXRATIO 1000.0

/* ------------------------------------------------------------
 * Constructors and destructors
 * ------------------------------------------------------------ */

pclustergeometry
new_clustergeometry(uint dim, uint nidx)
{
  pclustergeometry cg;
  real *buf;
  uint i;

  cg = (pclustergeometry) allocmem(sizeof(clustergeometry));
  cg->dim = dim;
  cg->nidx = nidx;

  cg->x = (real **) allocmem((size_t) sizeof(real *) * nidx);
  cg->smin = (real **) allocmem((size_t) sizeof(real *) * nidx);
  cg->smax = (real **) allocmem((size_t) sizeof(real *) * nidx);
  cg->perm = (uint *) allocmem((size_t) sizeof(uint) * nidx * (dim+2));

  cg->buf = buf = allocreal(dim * (3 * nidx + 2));

  cg->pflag = cg->perm + nidx * dim;
  cg->pbuf = cg->pflag + nidx;

  cg->pmin = buf;
  buf += dim;
  cg->pmax = buf;
  buf += dim;

  for(i=0; i<nidx; i++) {
    cg->x[i] = buf;
    buf += dim;
    cg->smin[i] = buf;
    buf += dim;
    cg->smax[i] = buf;
    buf += dim;
  }

  return cg;
}

void
del_clustergeometry(pclustergeometry cg)
{
  freemem(cg->buf);
  freemem(cg->perm);
  freemem(cg->smax);
  freemem(cg->smin);
  freemem(cg->x);
  freemem(cg);
}

/* ------------------------------------------------------------
 * Internal: Set bounding box for leaf cluster
 * ------------------------------------------------------------ */

static void
leafbbox(pcclustergeometry cg, pcluster t,
	 pcreal bmin, pcreal bmax)
{
  uint dim = t->dim;
  uint i, j, jj;

  assert(cg->dim == dim);

  for(i=0; i<dim; i++) {
    t->bmin[i] = bmin[i];
    t->bmax[i] = bmax[i];
  }

  for(j=0; j<t->size; j++) {
    jj = t->idx[j];
    for(i=0; i<dim; i++) {
      if(cg->smin[jj][i] < t->bmin[i])
	t->bmin[i] = cg->smin[jj][i];
      if(t->bmax[i] < cg->smax[jj][i])
	t->bmax[i] = cg->smax[jj][i];
    }
  }
}

/* ------------------------------------------------------------
 * Internal: Set bounding box for non-leaf cluster
 * ------------------------------------------------------------ */

static void
nonleafbbox(pcluster t,
	    pcreal bmin, pcreal bmax)
{
  uint dim = t->dim;
  uint i, j;

  for(i=0; i<dim; i++) {
    t->bmin[i] = bmin[i];
    t->bmax[i] = bmax[i];
  }

  for(j=0; j<t->sons; j++)
    for(i=0; i<dim; i++) {
      if(t->son[j]->bmin[i] < t->bmin[i])
	t->bmin[i] = t->son[j]->bmin[i];
      if(t->son[j]->bmax[i] > t->bmax[i])
	t->bmax[i] = t->son[j]->bmax[i];
  }
}

/* ------------------------------------------------------------
 * Internal: Fix large aspect ratios
 * ------------------------------------------------------------ */

static void
fixbbox(pcluster t, real maxratio)
{
  uint dim = t->dim;
  real width, maxwidth, middle;
  uint i;

  /* Find maximal width */
  maxwidth = t->bmax[0] - t->bmin[0];
  for(i=1; i<dim; i++) {
    width = t->bmax[i] - t->bmin[i];
    if(width > maxwidth)
      maxwidth = width;
  }

  /* Enlarge bounding box if maxwidth/width greater than maxratio */
  for(i=0; i<dim; i++) {
    width = t->bmax[i] - t->bmin[i];
    if(maxratio * width < maxwidth) {
      middle = 0.5 * (t->bmax[i] + t->bmin[i]);
      t->bmin[i] = middle - 0.5 * maxwidth / maxratio;
      t->bmax[i] = middle + 0.5 * maxwidth / maxratio;
    }
  }
}

/* ------------------------------------------------------------
 * Internal: Find bounding box for points
 * ------------------------------------------------------------ */

static void
findpbox(uint size, uint dim,
	 uint *idx, pcreal *x,
	 preal pmin, preal pmax)
{
  uint i, j;

  for(j=0; j<dim; j++)
    pmin[j] = pmax[j] = x[idx[0]][j];
  for(i=1; i<size; i++)
    for(j=0; j<dim; j++)
      if(pmin[j] > x[idx[i]][j])
	pmin[j] = x[idx[i]][j];
      else if(pmax[j] < x[idx[i]][j])
	pmax[j] = x[idx[i]][j];
}

/* ------------------------------------------------------------
 * Adaptive geometric cluster tree
 * ------------------------------------------------------------ */

pcluster
splitgeometric_clustergeometry(pclustergeometry cg,
			       int levels, uint resolution,
			       uint size, uint *idx)
{
  pcluster t;
  pcreal *x = (pcreal *) cg->x;
  preal pmin = cg->pmin;
  preal pmax = cg->pmax;
  real maxdiam, middle;
  uint dim = cg->dim;
  uint i, j, k, l;

  /* Find minimal bounding box for characteristic points */
  findpbox(size, dim, idx, x, pmin, pmax);
  
  if(levels <= 0 || size <= resolution) {
    /* Leaf cluster */
    t = new_cluster(size, idx, 0, dim);

    leafbbox(cg, t, pmin, pmax);
  }
  else {
    /* Find coordinate of maximal diameter */
    maxdiam = pmax[0] - pmin[0];
    k = 0;
    for(j=1; j<dim; j++)
      if(maxdiam < pmax[j] - pmin[j]) {
	maxdiam = pmax[j] - pmin[j];
	k = j;
      }
    
    /* Split index set */
    middle = 0.5 * (pmax[k] + pmin[k]);
    i = 0;
    j = size - 1;
    while(i < j) {
      while(i < size && x[idx[i]][k] <= middle)
	i++;
      
      while(0 < j && x[idx[j]][k] > middle)
	j--;
      
      if(i < j) {
	l = idx[i];
	idx[i] = idx[j];
	idx[j] = l;
      }
    }

    /* Father cluster */
    t = new_cluster(size, idx, 2, dim);

    /* Treat son clusters */
    t->son[0] = splitgeometric_clustergeometry(cg, levels-1, resolution,
					       i, idx);
    t->son[1] = splitgeometric_clustergeometry(cg, levels-1, resolution,
					       size-i, idx+i);

    nonleafbbox(t, pmin, pmax);
  }

  /* Fix excessively large aspect ratios */
  fixbbox(t, BBOX_MAXRATIO);

  /* Finish cluster */
  update_cluster(t);
  
  return t;
}

pcluster
buildgeometric_clustergeometry(pclustergeometry cg,
			       uint levels, uint resolution,
			       uint size, uint *idx)
{
  assert(resolution > 0);
  assert(size > 0);
  assert(idx != 0);

  return splitgeometric_clustergeometry(cg, levels, resolution,
					size, idx);
}

/* ------------------------------------------------------------
 * Regular geometric cluster tree
 * ------------------------------------------------------------ */

pcluster
splitregular_clustergeometry(pclustergeometry cg,
			     int levels, uint resolution, uint split,
			     uint size, uint *idx)
{
  pcluster t;
  pcreal *x = (pcreal *) cg->x;
  preal pmin = cg->pmin;
  preal pmax = cg->pmax;
  real original, middle;
  uint dim = cg->dim;
  uint i, j, l;

  assert(split < dim);

  if(levels <= 0 || size <= resolution) {
    /* Leaf cluster */
    t = new_cluster(size, idx, 0, dim);

    leafbbox(cg, t, pmin, pmax);
  }
  else {
    /* Split index set */
    middle = 0.5 * (pmax[split] + pmin[split]);
    i = 0;
    j = size - 1;
    while(i < j) {
      while(i < size && x[idx[i]][split] <= middle)
	i++;
      
      while(0 < j && x[idx[j]][split] > middle)
	j--;
      
      if(i < j) {
	l = idx[i];
	idx[i] = idx[j];
	idx[j] = l;
      }
    }
    
    /* Treat son clusters */
    if(i == size) {
      /* Only a "left" son */
      t = new_cluster(size, idx, 1, dim);

      original = pmax[split];
      pmax[split] = middle;
      t->son[0] = splitregular_clustergeometry(cg, levels-1, resolution,
					       (split+1)%dim, 
					       i, idx);
      pmax[split] = original;
    }
    else if(i == 0) {
      /* Only a "right" son */
      t = new_cluster(size, idx, 1, dim);

      original = pmin[split];
      pmin[split] = middle;
      t->son[0] = splitregular_clustergeometry(cg, levels-1, resolution,
					       (split+1)%dim,
					       size-i, idx);
      pmin[split] = original;
    }
    else {
      /* Two sons */
      t = new_cluster(size, idx, 2, dim);
      
      original = pmax[split];
      pmax[split] = middle;
      t->son[0] = splitregular_clustergeometry(cg, levels-1, resolution,
					       (split+1)%dim,
					       i, idx);
      pmax[split] = original;

      original = pmin[split];
      pmin[split] = middle;
      t->son[1] = splitregular_clustergeometry(cg, levels-1, resolution,
					       (split+1)%dim,
					       size-i, idx+i);
      pmin[split] = original;
    }

    nonleafbbox(t, pmin, pmax);
  }

  /* Fix excessively large aspect ratios */
  fixbbox(t, BBOX_MAXRATIO);

  /* Finish cluster */
  update_cluster(t);
  
  return t;
}

pcluster
buildregular_clustergeometry(pclustergeometry cg,
			     uint levels, uint resolution,
			     uint size, uint *idx)
{
  assert(resolution > 0);
  assert(size > 0);
  assert(idx != 0);

  findpbox(size, cg->dim, idx, (pcreal *) cg->x, cg->pmin, cg->pmax);

  return splitregular_clustergeometry(cg, levels, resolution, 0,
				      size, idx);
}

/* ------------------------------------------------------------
 * Cardinality-balanced cluster tree
 * ------------------------------------------------------------ */

pcluster
splitcardinality_clustergeometry(pclustergeometry cg,
				 int levels, uint resolution,
				 uint size, uint *idx,
				 uint *perm)
{
  pcluster t;
  pcreal *x = (pcreal *) cg->x;
  preal pmin = cg->pmin;
  preal pmax = cg->pmax;
  uint *pflag = cg->pflag;
  uint *pbuf = cg->pbuf;
  real maxdiam;
  uint dim = cg->dim;
  uint i, j, k, l;

  /* Find minimal bounding box for characteristic points */
  findpbox(size, dim, idx, x, pmin, pmax);
  
  if(levels <= 0 || size <= resolution) {
    /* Leaf cluster */
    t = new_cluster(size, idx, 0, dim);

    leafbbox(cg, t, pmin, pmax);
  }
  else {
    /* Find coordinate of maximal diameter */
    maxdiam = pmax[0] - pmin[0];
    k = 0;
    for(j=1; j<dim; j++)
      if(maxdiam < pmax[j] - pmin[j]) {
	maxdiam = pmax[j] - pmin[j];
	k = j;
      }

    /* Use k-th ordered array */
    for(i=0; i<size; i++)
      idx[i] = perm[k+i*dim];
    
    /* Mark indices in left and right cluster */
    for(i=0; i<size/2; i++)
      pflag[perm[k+i*dim]] = 1;
    for(; i<size; i++)
      pflag[perm[k+i*dim]] = 0;

    /* Split sorted arrays for each coordinate to create
       sorted arrays for left and right clusters */
    for(j=0; j<dim; j++) {
      l = 0;

      for(i=0; i<size; i++)
	if(pflag[perm[j+i*dim]]) {
	  pbuf[l] = perm[j+i*dim];
	  l++;
	}
      assert(l == size/2);

      for(i=0; i<size; i++)
	if(!pflag[perm[j+i*dim]]) {
	  pbuf[l] = perm[j+i*dim];
	  l++;
	}
      assert(l == size);

      for(i=0; i<size; i++)
	perm[j+i*dim] = pbuf[i];
    }

    /* Father cluster */
    t = new_cluster(size, idx, 2, dim);

    /* Treat son clusters */
    t->son[0] = splitcardinality_clustergeometry(cg, levels-1, resolution,
						 size/2, idx,
						 perm);
    t->son[1] = splitcardinality_clustergeometry(cg, levels-1, resolution,
						 size-size/2, idx+size/2,
						 perm+(size/2)*dim);

    nonleafbbox(t, pmin, pmax);
  }

  /* Fix excessively large aspect ratios */
  fixbbox(t, BBOX_MAXRATIO);

  /* Finish cluster */
  update_cluster(t);
  
  return t;
}

struct _coordsort {
  pcreal *x;
  uint *idx;
  uint k;
};

static bool coordsort_leq(uint i, uint j, void *sd)
{
  struct _coordsort *cs = (struct _coordsort *) sd;

  return (cs->x[cs->idx[i]][cs->k] <= cs->x[cs->idx[j]][cs->k]);
}

static void coordsort_swap(uint i, uint j, void *sd)
{
  struct _coordsort *cs = (struct _coordsort *) sd;
  uint b;

  b = cs->idx[i];
  cs->idx[i] = cs->idx[j];
  cs->idx[j] = b;
}

pcluster
buildcardinality_clustergeometry(pclustergeometry cg,
				 uint levels, uint resolution,
				 uint size, uint *idx)
{
  struct _coordsort cs;
  pcreal *x = (pcreal *) cg->x;
  uint *perm = cg->perm;
  uint *pbuf = cg->pbuf;
  uint dim = cg->dim;
  uint i, k;

  assert(resolution > 0);
  assert(size > 0);
  assert(idx != 0);

  /* Sort coordinates in all directions */
  cs.x = (pcreal *) x;
  for(k=0; k<dim; k++) {
    for(i=0; i<size; i++)
      pbuf[i] = idx[i];
    cs.idx = pbuf;
    cs.k = k;

    heapsort(size, coordsort_leq, coordsort_swap, &cs);

    for(i=0; i<size; i++)
      perm[k+i*dim] = pbuf[i];

#ifndef NDEBUG
    for(i=0; i<size-1; i++)
      assert(x[perm[k+i*dim]][k] <= x[perm[k+(i+1)*dim]][k]);
#endif
  }

  return splitcardinality_clustergeometry(cg, levels, resolution,
					  size, idx, perm);
}

/* ------------------------------------------------------------
 * Cardinality/geometric cluster tree
 * ------------------------------------------------------------ */

pcluster
splitmixed_clustergeometry(pclustergeometry cg,
			   uint cardepth, int levels, uint resolution,
			   uint size, uint *idx,
			   uint *perm)
{
  pcluster t;
  pcreal *x = (pcreal *) cg->x;
  preal pmin = cg->pmin;
  preal pmax = cg->pmax;
  uint *pflag = cg->pflag;
  uint *pbuf = cg->pbuf;
  real maxdiam;
  uint dim = cg->dim;
  uint i, j, k, l;

  /* Find minimal bounding box for characteristic points */
  findpbox(size, dim, idx, x, pmin, pmax);
  
  if(levels <= 0 || size <= resolution) {
    /* Leaf cluster */
    t = new_cluster(size, idx, 0, dim);

    leafbbox(cg, t, pmin, pmax);
  }
  else {
    /* Find coordinate of maximal diameter */
    maxdiam = pmax[0] - pmin[0];
    k = 0;
    for(j=1; j<dim; j++)
      if(maxdiam < pmax[j] - pmin[j]) {
	maxdiam = pmax[j] - pmin[j];
	k = j;
      }

    /* Use k-th ordered array */
    for(i=0; i<size; i++)
      idx[i] = perm[k+i*dim];
    
    /* Mark indices in left and right cluster */
    for(i=0; i<size/2; i++)
      pflag[perm[k+i*dim]] = 1;
    for(; i<size; i++)
      pflag[perm[k+i*dim]] = 0;

    /* Split sorted arrays for each coordinate to create
       sorted arrays for left and right clusters */
    for(j=0; j<dim; j++) {
      l = 0;

      for(i=0; i<size; i++)
	if(pflag[perm[j+i*dim]]) {
	  pbuf[l] = perm[j+i*dim];
	  l++;
	}
      assert(l == size/2);

      for(i=0; i<size; i++)
	if(!pflag[perm[j+i*dim]]) {
	  pbuf[l] = perm[j+i*dim];
	  l++;
	}
      assert(l == size);

      for(i=0; i<size; i++)
	perm[j+i*dim] = pbuf[i];
    }

    /* Father cluster */
    t = new_cluster(size, idx, 2, dim);

    /* Treat son clusters */
    if(cardepth > 0) {
      t->son[0] = splitmixed_clustergeometry(cg, cardepth-1,
					     levels-1, resolution,
					     size/2, idx,
					     perm);
      t->son[1] = splitmixed_clustergeometry(cg, cardepth-1,
					     levels-1, resolution,
					     size-size/2, idx+size/2,
					     perm+(size/2)*dim);
    }
    else {
      t->son[0] = splitgeometric_clustergeometry(cg, levels-1, resolution,
						 size/2, idx);
      t->son[1] = splitgeometric_clustergeometry(cg, levels-1, resolution,
						 size-size/2, idx+size/2);
    }

    nonleafbbox(t, pmin, pmax);
  }
  
  /* Fix excessively large aspect ratios */
  fixbbox(t, BBOX_MAXRATIO);

  /* Finish cluster */
  update_cluster(t);
  
  return t;
}

pcluster
buildmixed_clustergeometry(pclustergeometry cg, uint cardepth,
			   uint levels, uint resolution,
			   uint size, uint *idx)
{
  struct _coordsort cs;
  pcreal *x = (pcreal *) cg->x;
  uint *perm = cg->perm;
  uint *pbuf = cg->pbuf;
  uint dim = cg->dim;
  uint i, k;

  assert(resolution > 0);
  assert(size > 0);
  assert(idx != 0);

  /* Sort coordinates in all directions */
  cs.x = (pcreal *) x;
  for(k=0; k<dim; k++) {
    for(i=0; i<size; i++)
      pbuf[i] = idx[i];
    cs.idx = pbuf;
    cs.k = k;

    heapsort(size, coordsort_leq, coordsort_swap, &cs);

    for(i=0; i<size; i++)
      perm[k+i*dim] = pbuf[i];

#ifndef NDEBUG
    for(i=0; i<size-1; i++)
      assert(x[perm[k+i*dim]][k] <= x[perm[k+(i+1)*dim]][k]);
#endif
  }

  return splitmixed_clustergeometry(cg, cardepth, levels, resolution,
				    size, idx, perm);
}
