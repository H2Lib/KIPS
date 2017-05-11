
/* ------------------------------------------------------------
 * This is the file "clusterbasis.c" of the KIPS package.
 * All rights reserved, Steffen Boerm 2009
 * ------------------------------------------------------------ */

#include "clusterbasis.h"

#include "basic.h"

static uint active_clusterbasis = 0;

/* ------------------------------------------------------------
 * Constructors and destructors
 * ------------------------------------------------------------ */

static void
init_raw_clusterbasis(pclusterbasis cb, pccluster t)
{
  assert(cb != NULL);

  cb->t = t;
  cb->k = 0;
  cb->ktree = 0;

  cb->sons = 0;
  cb->son = 0;

#ifdef USE_OPENMP
#pragma omp atomic
#endif
  active_clusterbasis++;
}

pclusterbasis
init_clusterbasis(pclusterbasis cb, pccluster t)
{
  uint      i, sons;

  (void) init_leaf_clusterbasis(cb, t);

  sons = t->sons;
  if (sons > 0) {
    cb->sons = sons;
    cb->son =
      (pclusterbasis *) allocmem((size_t) sizeof(pclusterbasis) * sons);
    for (i = 0; i < sons; i++)
      cb->son[i] = 0;
  }

  return cb;
}

pclusterbasis
init_leaf_clusterbasis(pclusterbasis cb, pccluster t)
{
  assert(cb != NULL);

  init_raw_clusterbasis(cb, t);

  init_amatrix(&cb->V, 0, 0);
  init_amatrix(&cb->E, 0, 0);

  return cb;
}

void
uninit_clusterbasis(pclusterbasis cb)
{
  uint      i;

  if (cb->sons > 0) {
    for (i = 0; i < cb->sons; i++)
      del_clusterbasis(cb->son[i]);
    freemem(cb->son);
  }

  uninit_amatrix(&cb->V);
  uninit_amatrix(&cb->E);

  assert(active_clusterbasis > 0);

#ifdef USE_OPENMP
#pragma omp atomic
#endif
  active_clusterbasis--;
}

pclusterbasis
new_clusterbasis(pccluster t)
{
  pclusterbasis cb;

  cb = allocmem(sizeof(clusterbasis));

  init_clusterbasis(cb, t);

  return cb;
}

pclusterbasis
new_leaf_clusterbasis(pccluster t)
{
  pclusterbasis cb;

  cb = allocmem(sizeof(clusterbasis));

  init_leaf_clusterbasis(cb, t);

  return cb;
}

void
del_clusterbasis(pclusterbasis cb)
{
  uninit_clusterbasis(cb);

  freemem(cb);
}

/* ------------------------------------------------------------
 * Low-level management
 * ------------------------------------------------------------ */

void
update_clusterbasis(pclusterbasis cb)
{
  uint      stree;
  uint      i;

  stree = 0;

  if (cb->sons == 0) {
    stree += cb->t->size;
  }
  else {
    for (i = 0; i < cb->sons; i++)
      stree += cb->son[i]->ktree;
  }

  cb->ktree = stree + cb->k;
}

void
setrank_clusterbasis(uint k, pclusterbasis cb)
{
  uint      i;

  if (cb->sons > 0) {
    for (i = 0; i < cb->sons; i++)
      resize_amatrix(&cb->son[i]->E, cb->son[i]->k, k);
  }
  else
    resize_amatrix(&cb->V, cb->t->size, k);

  cb->k = k;

  update_clusterbasis(cb);
}

/* ------------------------------------------------------------
 * Build clusterbasis based on cluster
 * ------------------------------------------------------------ */

pclusterbasis
build_fromcluster_clusterbasis(pccluster t)
{
  pclusterbasis cb, cb1;
  uint      i;

  cb = new_clusterbasis(t);

  if (cb->sons > 0) {
    for (i = 0; i < cb->sons; i++) {
      cb1 = build_fromcluster_clusterbasis(t->son[i]);
      cb->son[i] = cb1;
    }
  }

  update_clusterbasis(cb);

  return cb;
}

/* ------------------------------------------------------------
 * Initialize clusterbasis using callback functions
 * ------------------------------------------------------------ */

void
fill_clusterbasis(void (*buildV)(void *data, pccluster t, pamatrix V),
		  void (*buildE)(void *data, pccluster sc, pccluster fc, pamatrix E),
		  void *data,
		  pclusterbasis cb)
{
  uint i;

  if(cb->sons > 0) {
    for(i=0; i<cb->sons; i++) {
      buildE(data, cb->son[i]->t, cb->t, &cb->son[i]->E);

      fill_clusterbasis(buildV, buildE, data, cb->son[i]);

      assert(i == 0 || cb->son[i]->E.cols == cb->son[0]->E.cols);
    }
    cb->k = cb->son[0]->E.cols;
  }
  else {
    buildV(data, cb->t, &cb->V);

    cb->k = cb->V.cols;
  }

  update_clusterbasis(cb);
}

/* ------------------------------------------------------------
 * Statistics
 * ------------------------------------------------------------ */

uint
getactives_clusterbasis()
{
  return active_clusterbasis;
}

size_t
getsize_clusterbasis(pcclusterbasis cb)
{
  size_t    sz;
  uint      i;

  sz = (size_t) sizeof(clusterbasis);
  sz += getsize_heap_amatrix(&cb->V);
  sz += getsize_heap_amatrix(&cb->E);

  if (cb->sons > 0) {
    sz += (size_t) sizeof(pclusterbasis) * cb->sons;

    for (i = 0; i < cb->sons; i++)
      sz += getsize_clusterbasis(cb->son[i]);
  }

  return sz;
}

/* ------------------------------------------------------------
 * Forward and backward transformation
 * ------------------------------------------------------------ */

void
forward_clusterbasis(pcclusterbasis cb, pcavector x, pavector xt)
{
  avector   loc1, loc2;
  pavector  xt1, xc, xp;
  uint      i, xtoff;

  assert(xt->size == cb->ktree);

  /* This part of xt contains the coefficients for the current cluster */
  xc = init_sub_avector(&loc1, xt, cb->k, 0);
  clear_avector(xc);

  if (cb->sons > 0) {
    xtoff = cb->k;
    for (i = 0; i < cb->sons; i++) {
      /* This part corresponds to the subtree rooted in the i-th son */
      xt1 = init_sub_avector(&loc2, xt, cb->son[i]->ktree, xtoff);

      /* Compute coefficients in the subtree */
      forward_clusterbasis(cb->son[i], x, xt1);

      uninit_avector(xt1);

      /* This part corresponds to the i-th son */
      xt1 = init_sub_avector(&loc2, xt, cb->son[i]->k, xtoff);

      /* Multiply by transfer matrix */
      addevaltrans_amatrix(1.0, &cb->son[i]->E, xt1, xc);

      uninit_avector(xt1);

      xtoff += cb->son[i]->ktree;
    }
    assert(xtoff == cb->ktree);
  }
  else {
    /* Permuted entries of x */
    xp = init_sub_avector(&loc2, xt, cb->t->size, cb->k);

    /* Find and copy entries */
    for (i = 0; i < cb->t->size; i++)
      xp->v[i] = x->v[cb->t->idx[i]];

    /* Multiply by leaf matrix */
    addevaltrans_amatrix(1.0, &cb->V, xp, xc);

    uninit_avector(xp);
  }

  uninit_avector(xc);
}

void
backward_clusterbasis(pcclusterbasis cb, pavector yt, pavector y)
{
  avector   loc1, loc2;
  pavector  yt1, yc, yp;
  uint      i, ytoff;

  assert(yt->size == cb->ktree);

  /* This part of yt contains the coefficients for the current cluster */
  yc = init_sub_avector(&loc1, yt, cb->k, 0);

  if (cb->sons > 0) {
    ytoff = cb->k;
    for (i = 0; i < cb->sons; i++) {
      /* This part corresponds to the i-th son */
      yt1 = init_sub_avector(&loc2, yt, cb->son[i]->k, ytoff);

      /* Multiply by transfer matrix */
      addeval_amatrix(1.0, &cb->son[i]->E, yc, yt1);

      uninit_avector(yt1);

      /* This part corresponds to the subtree rooted in the i-th son */
      yt1 = init_sub_avector(&loc2, yt, cb->son[i]->ktree, ytoff);

      /* Treat coefficients in the subtree */
      backward_clusterbasis(cb->son[i], yt1, y);

      uninit_avector(yt1);

      ytoff += cb->son[i]->ktree;
    }
    assert(ytoff == cb->ktree);
  }
  else {
    /* Permuted entries of x */
    yp = init_sub_avector(&loc2, yt, cb->t->size, cb->k);

    /* Multiply by leaf matrix */
    addeval_amatrix(1.0, &cb->V, yc, yp);

    /* Find and copy entries */
    for (i = 0; i < cb->t->size; i++)
      y->v[cb->t->idx[i]] += yp->v[i];

    uninit_avector(yp);
  }

  uninit_avector(yc);
}
