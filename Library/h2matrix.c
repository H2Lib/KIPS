
/* ------------------------------------------------------------
 * This is the file "h2matrix.c" of the KIPS package.
 * All rights reserved, Steffen Boerm 2009
 * ------------------------------------------------------------ */

#include <stdio.h>

#include "h2matrix.h"

#include "basic.h"

/* ------------------------------------------------------------
 * Constructors and destructors
 * ------------------------------------------------------------ */

ph2matrix
new_h2matrix(pclusterbasis rb, pclusterbasis cb)
{
  ph2matrix h2;

  h2 = allocmem(sizeof(h2matrix));

  h2->rb = rb;
  h2->cb = cb;

  h2->u = NULL;
  h2->f = NULL;

  h2->son = NULL;
  h2->rsons = 0;
  h2->csons = 0;

  h2->desc = 0;

  return h2;
}

ph2matrix
new_uniform_h2matrix(pclusterbasis rb, pclusterbasis cb)
{
  ph2matrix h2;

  h2 = new_h2matrix(rb, cb);

  h2->u = new_uniform(rb, cb);

  h2->desc = 1;

  return h2;
}

ph2matrix
new_full_h2matrix(pclusterbasis rb, pclusterbasis cb)
{
  ph2matrix h2;

  h2 = new_h2matrix(rb, cb);

  h2->f = new_amatrix(rb->t->nidx, cb->t->nidx);

  h2->desc = 1;

  return h2;
}

ph2matrix
new_super_h2matrix(pclusterbasis rb, pclusterbasis cb, uint rsons, uint csons)
{
  ph2matrix h2;
  uint      i, j;

  h2 = new_h2matrix(rb, cb);

  h2->rsons = rsons;
  h2->csons = csons;

  h2->son =
    (ph2matrix *) allocmem((size_t) sizeof(ph2matrix) * rsons * csons);
  for (j = 0; j < csons; j++)
    for (i = 0; i < rsons; i++)
      h2->son[i + j * rsons] = NULL;

  return h2;
}

ph2matrix
new_zero_h2matrix(pclusterbasis rb, pclusterbasis cb)
{
  ph2matrix h2;

  h2 = new_h2matrix(rb, cb);

  h2->u = 0;
  h2->f = 0;
  h2->son = 0;
  h2->desc = 1;

  return h2;
}

ph2matrix
clonestructure_h2matrix(pch2matrix G, pclusterbasis rb, pclusterbasis cb)
{
  ph2matrix H, H1;
  pch2matrix G1;
  pclusterbasis rb1, cb1;
  uint      rsons, csons;
  uint      i, j;

  H = NULL;

  if (G->son) {
    rsons = G->rsons;
    csons = G->csons;

    H = new_super_h2matrix(rb, cb, rsons, csons);

    for (j = 0; j < csons; j++)
      for (i = 0; i < rsons; i++) {
	G1 = G->son[i + j * rsons];

	rb1 = rb;
	if (G1->rb->t != G->rb->t) {
	  assert(rb->sons == rsons);
	  rb1 = rb->son[i];
	}

	cb1 = cb;
	if (G1->cb->t != G->cb->t) {
	  assert(cb->sons == csons);
	  cb1 = cb->son[j];
	}

	H1 = clonestructure_h2matrix(G1, rb1, cb1);

	H->son[i + j * rsons] = H1;
      }

  }
  else if (G->u) {
    H = new_uniform_h2matrix(rb, cb);
    clear_amatrix(&H->u->S);
  }
  else if (G->f) {
    H = new_full_h2matrix(rb, cb);
    clear_amatrix(H->f);
  }
  else
    H = new_zero_h2matrix(rb, cb);

  update_h2matrix(H);

  return H;
}

ph2matrix
clone_h2matrix(pch2matrix G, pclusterbasis rb, pclusterbasis cb)
{
  ph2matrix H, H1;
  pch2matrix G1;
  pclusterbasis rb1, cb1;
  uint      rsons, csons;
  uint      i, j;

  H = NULL;

  if (G->son) {
    rsons = G->rsons;
    csons = G->csons;

    H = new_super_h2matrix(rb, cb, rsons, csons);

    for (j = 0; j < csons; j++)
      for (i = 0; i < rsons; i++) {
	G1 = G->son[i + j * rsons];

	rb1 = rb;
	if (G1->rb->t != G->rb->t) {
	  assert(rb->sons == rsons);
	  rb1 = rb->son[i];
	}

	cb1 = cb;
	if (G1->cb->t != G->cb->t) {
	  assert(cb->sons == csons);
	  cb1 = cb->son[j];
	}

	H1 = clone_h2matrix(G1, rb1, cb1);

	H->son[i + j * rsons] = H1;
      }

  }
  else if (G->u) {
    H = new_uniform_h2matrix(rb, cb);
    copy_uniform(false, G->u, H->u);
  }
  else if (G->f) {
    H = new_full_h2matrix(rb, cb);
    copy_amatrix(false, G->f, H->f);
  }
  else
    H = new_zero_h2matrix(rb, cb);

  update_h2matrix(H);

  return H;
}

void
update_h2matrix(ph2matrix h2)
{
  uint      desc;
  uint      rsons, csons;
  uint      i, j;

  desc = 1;

  if (h2->son) {
    rsons = h2->rsons;
    csons = h2->csons;

    for (j = 0; j < csons; j++)
      for (i = 0; i < rsons; i++)
	desc += h2->son[i + j * rsons]->desc;
  }

  h2->desc = desc;
}

void
del_h2matrix(ph2matrix h2)
{
  uint      rsons = h2->rsons;
  uint      csons = h2->csons;
  uint      i, j;

  if (h2->son) {
    for (j = 0; j < csons; j++)
      for (i = 0; i < rsons; i++)
	del_h2matrix(h2->son[i + j * rsons]);
    freemem(h2->son);
  }

  if (h2->f)
    del_amatrix(h2->f);

  if (h2->u)
    del_uniform(h2->u);

  h2->cb = 0;
  h2->rb = 0;

  freemem(h2);
}

/* ------------------------------------------------------------
 * Statistics
 * ------------------------------------------------------------ */

size_t
getsize_h2matrix(pch2matrix h2)
{
  size_t    sz;
  uint      rsons = h2->rsons;
  uint      csons = h2->csons;
  uint      i, j;

  sz = (size_t) sizeof(h2matrix);

  if (h2->u)
    sz += getsize_uniform(h2->u);

  if (h2->f)
    sz += getsize_amatrix(h2->f);

  for (j = 0; j < csons; j++)
    for (i = 0; i < rsons; i++)
      sz += getsize_h2matrix(h2->son[i + j * rsons]);

  return sz;
}

size_t
gettotalsize_h2matrix(pch2matrix h2)
{
  size_t    sz;

  sz = getsize_h2matrix(h2);
  sz += getsize_clusterbasis(h2->rb);
  if (h2->rb != h2->cb)
    sz += getsize_clusterbasis(h2->cb);

  return sz;
}

size_t
getnearsize_h2matrix(pch2matrix h2)
{
  size_t    sz;
  uint      rsons = h2->rsons;
  uint      csons = h2->csons;
  uint      i, j;

  sz = 0;

  if (h2->f)
    sz += getsize_amatrix(h2->f);

  for (j = 0; j < csons; j++)
    for (i = 0; i < rsons; i++)
      sz += getnearsize_h2matrix(h2->son[i + j * rsons]);

  return sz;
}

size_t
getfarsize_h2matrix(pch2matrix h2)
{
  size_t    sz;
  uint      rsons = h2->rsons;
  uint      csons = h2->csons;
  uint      i, j;

  sz = 0;

  if (h2->u)
    sz += getsize_uniform(h2->u);

  for (j = 0; j < csons; j++)
    for (i = 0; i < rsons; i++)
      sz += getfarsize_h2matrix(h2->son[i + j * rsons]);

  return sz;
}

/* ------------------------------------------------------------
 * Simple utility functions
 * ------------------------------------------------------------ */

void
clear_h2matrix(ph2matrix h2)
{
  uint      rsons, csons;
  uint      i, j;

  if (h2->son) {
    rsons = h2->rsons;
    csons = h2->csons;

    for (j = 0; j < csons; j++)
      for (i = 0; i < rsons; i++)
	clear_h2matrix(h2->son[i + j * rsons]);
  }
  else if (h2->u)
    clear_amatrix(&h2->u->S);
  else if (h2->f)
    clear_amatrix(h2->f);
}

void
scale_h2matrix(field alpha, ph2matrix h2)
{
  uint      rsons = h2->rsons;
  uint      csons = h2->csons;
  uint      i, j;

  if (h2->son) {
    for (j = 0; j < csons; j++)
      for (i = 0; i < rsons; i++)
	scale_h2matrix(alpha, h2->son[i + j * rsons]);
  }
  else if (h2->u)
    scale_uniform(alpha, h2->u);
  else if (h2->f)
    scale_amatrix(alpha, h2->f);
}

void
random_h2matrix(ph2matrix h2)
{
  uint      rsons = h2->rsons;
  uint      csons = h2->csons;
  uint      i, j;

  if (h2->son) {
    for (j = 0; j < csons; j++)
      for (i = 0; i < rsons; i++)
	random_h2matrix(h2->son[i + j * rsons]);
  }
  else if (h2->u)
    random_uniform(h2->u);
  else if (h2->f)
    random_amatrix(h2->f);
}

/* ------------------------------------------------------------
 * Build H^2-matrix based on block tree
 * ------------------------------------------------------------ */

ph2matrix
build_fromblock_h2matrix(pcblock b, pclusterbasis rb, pclusterbasis cb)
{
  ph2matrix h, h1;
  pcblock   b1;
  pclusterbasis rb1, cb1;
  uint      rsons, csons;
  uint      i, j;

  h = NULL;

  if (b->son) {
    rsons = b->rsons;
    csons = b->csons;

    h = new_super_h2matrix(rb, cb, rsons, csons);

    for (j = 0; j < csons; j++)
      for (i = 0; i < rsons; i++) {
	b1 = b->son[i + j * rsons];

	rb1 = rb;
	if (b1->rc != b->rc) {
	  assert(rb->sons == rsons);
	  rb1 = rb->son[i];
	}

	cb1 = cb;
	if (b1->cc != b->cc) {
	  assert(cb->sons == csons);
	  cb1 = cb->son[j];
	}

	h1 = build_fromblock_h2matrix(b1, rb1, cb1);

	h->son[i + j * rsons] = h1;
      }
  }
  else if (b->adm)
    h = new_uniform_h2matrix(rb, cb);
  else
    h = new_full_h2matrix(rb, cb);

  update_h2matrix(h);

  return h;
}

/* ------------------------------------------------------------
 * Initialize H^2-matrix using callback functions
 * ------------------------------------------------------------ */

void
fill_h2matrix(void (*buildN)
	      (void *data, const uint * ridx, const uint * cidx, pamatrix N),
	      void (*buildS) (void *data, pcspatialcluster rc, pcspatialcluster cc,
			      pamatrix S), void *data, ph2matrix G)
{
  uint      rsons, csons;
  uint      i, j;

  if (G->son) {
    rsons = G->rsons;
    csons = G->csons;

    for (j = 0; j < csons; j++)
      for (i = 0; i < rsons; i++)
	fill_h2matrix(buildN, buildS, data, G->son[i + j * rsons]);
  }
  else if (G->u) {
    buildS(data, G->rb->t, G->cb->t, &G->u->S);

    assert(G->u->S.rows == G->rb->k);
    assert(G->u->S.cols == G->cb->k);
  }
  else if (G->f) {
    buildN(data, G->rb->t->idx, G->cb->t->idx, G->f);

    assert(G->f->rows == G->rb->t->nidx);
    assert(G->f->cols == G->cb->t->nidx);
  }
}

/* ------------------------------------------------------------
 * Matrix-vector multiplication
 * ------------------------------------------------------------ */

void
mvm_h2matrix(field alpha, bool h2trans, pch2matrix h2,
	     pcavector x, pavector y)
{
  if (h2trans)
    addevaltrans_h2matrix(alpha, h2, x, y);
  else
    addeval_h2matrix(alpha, h2, x, y);
}

void
fastaddeval_h2matrix(field alpha, pch2matrix h2, pcavector xt, pavector yt)
{
  avector   loc1, loc2;
  pavector  xp, yp, xt1, yt1;
  pcclusterbasis rb = h2->rb;
  pcclusterbasis cb = h2->cb;
  uint      rsons = h2->rsons;
  uint      csons = h2->csons;
  uint      xtoff, ytoff;
  uint      i, j;

  assert(xt->size == h2->cb->ktree);
  assert(yt->size == h2->rb->ktree);

  if (h2->u) {
    xt1 = init_sub_avector(&loc1, (pavector) xt, cb->k, 0);
    yt1 = init_sub_avector(&loc2, yt, rb->k, 0);

    addeval_amatrix(alpha, &h2->u->S, xt1, yt1);

    uninit_avector(yt1);
    uninit_avector(xt1);
  }
  else if (h2->f) {
    xp = init_sub_avector(&loc1, (pavector) xt, cb->t->nidx, cb->k);
    yp = init_sub_avector(&loc2, yt, rb->t->nidx, rb->k);

    addeval_amatrix(alpha, h2->f, xp, yp);

    uninit_avector(yp);
    uninit_avector(xp);
  }
  else if (h2->son) {
    xtoff = cb->k;
    for (j = 0; j < csons; j++) {
      assert(csons == 1 || cb->sons > 0);
      xt1 =
	(cb->sons > 0 ?
	 init_sub_avector(&loc1, (pavector) xt, cb->son[j]->ktree, xtoff) :
	 init_sub_avector(&loc1, (pavector) xt, cb->ktree, 0));

      ytoff = rb->k;
      for (i = 0; i < rsons; i++) {
	assert(rsons == 1 || rb->sons > 0);
	yt1 = (rb->sons > 0 ?
	       init_sub_avector(&loc2, yt, rb->son[i]->ktree, ytoff) :
	       init_sub_avector(&loc2, yt, rb->ktree, 0));

	fastaddeval_h2matrix(alpha, h2->son[i + j * rsons], xt1, yt1);

	uninit_avector(yt1);

	ytoff += (rb->sons > 0 ? rb->son[i]->ktree : rb->t->nidx);
      }
      assert(ytoff == rb->ktree);

      uninit_avector(xt1);

      xtoff += (cb->sons > 0 ? cb->son[j]->ktree : cb->t->nidx);
    }
    assert(xtoff == cb->ktree);
  }
}

void
addeval_h2matrix(field alpha, pch2matrix h2, pcavector x, pavector y)
{
  pavector  xt, yt;

  xt = new_avector(h2->cb->ktree);
  yt = new_avector(h2->rb->ktree);

  clear_avector(yt);

  forward_clusterbasis(h2->cb, x, xt);

  fastaddeval_h2matrix(alpha, h2, xt, yt);

  backward_clusterbasis(h2->rb, yt, y);

  del_avector(yt);
  del_avector(xt);
}

void
fastaddevaltrans_h2matrix(field alpha, pch2matrix h2,
			  pcavector xt, pavector yt)
{
  avector   loc1, loc2;
  pavector  xp, yp, xt1, yt1;
  pcclusterbasis rb = h2->rb;
  pcclusterbasis cb = h2->cb;
  uint      rsons = h2->rsons;
  uint      csons = h2->csons;
  uint      xtoff, ytoff;
  uint      i, j;

  assert(xt->size == h2->rb->ktree);
  assert(yt->size == h2->cb->ktree);

  if (h2->u) {
    xt1 = init_sub_avector(&loc1, (pavector) xt, rb->k, 0);
    yt1 = init_sub_avector(&loc2, yt, cb->k, 0);

    addevaltrans_amatrix(alpha, &h2->u->S, xt1, yt1);

    uninit_avector(yt1);
    uninit_avector(xt1);
  }
  else if (h2->f) {
    xp = init_sub_avector(&loc1, (pavector) xt, rb->t->nidx, rb->k);
    yp = init_sub_avector(&loc2, yt, cb->t->nidx, cb->k);

    addevaltrans_amatrix(alpha, h2->f, xp, yp);

    uninit_avector(yp);
    uninit_avector(xp);
  }
  else if (h2->son) {
    ytoff = cb->k;
    for (j = 0; j < csons; j++) {
      yt1 =
	(cb->sons >
	 0 ? init_sub_avector(&loc2, yt, cb->son[j]->ktree,
			      ytoff) : init_sub_avector(&loc2, yt, cb->ktree,
							0));

      xtoff = rb->k;
      for (i = 0; i < rsons; i++) {
	xt1 = (rb->sons > 0 ?
	       init_sub_avector(&loc1, (pavector) xt, rb->son[i]->ktree,
				xtoff) : init_sub_avector(&loc1,
							  (pavector) xt,
							  rb->ktree, 0));

	fastaddevaltrans_h2matrix(alpha, h2->son[i + j * rsons], xt1, yt1);

	uninit_avector(xt1);

	xtoff += (rb->sons > 0 ? rb->son[i]->ktree : rb->t->nidx);
      }
      assert(xtoff == rb->ktree);

      uninit_avector(yt1);

      ytoff += (cb->sons > 0 ? cb->son[j]->ktree : cb->t->nidx);
    }
    assert(ytoff == cb->ktree);
  }
}

void
addevaltrans_h2matrix(field alpha, pch2matrix h2, pcavector x, pavector y)
{
  pavector  xt, yt;

  xt = new_avector(h2->rb->ktree);
  yt = new_avector(h2->cb->ktree);

  clear_avector(yt);

  forward_clusterbasis(h2->rb, x, xt);

  fastaddevaltrans_h2matrix(alpha, h2, xt, yt);

  backward_clusterbasis(h2->cb, yt, y);

  del_avector(yt);
  del_avector(xt);
}

static void
addevalsymm_offdiag(field alpha, pch2matrix h2, pavector xt,
		    pavector xta, pavector yt, pavector yta)
{
  avector   tmp1, tmp2, tmp3, tmp4;
  pavector  xp, yp;
  pavector  xt1, xta1, yt1, yta1;
  pcclusterbasis rb = h2->rb;
  pcclusterbasis cb = h2->cb;
  uint      rsons, csons;
  uint      xtoff, ytoff;
  uint      i, j;

  assert(xt->size == h2->cb->ktree);
  assert(xta->size == h2->rb->ktree);
  assert(yt->size == h2->rb->ktree);
  assert(yta->size == h2->cb->ktree);

  if (h2->f) {
    xp = init_sub_avector(&tmp1, xt, cb->t->nidx, cb->k);
    yp = init_sub_avector(&tmp2, yt, rb->t->nidx, rb->k);

    addeval_amatrix(alpha, h2->f, xp, yp);

    uninit_avector(yp);
    uninit_avector(xp);

    xp = init_sub_avector(&tmp1, xta, rb->t->nidx, rb->k);
    yp = init_sub_avector(&tmp2, yta, cb->t->nidx, cb->k);

    addevaltrans_amatrix(alpha, h2->f, xp, yp);

    uninit_avector(yp);
    uninit_avector(xp);
  }
  else if (h2->u) {
    addeval_amatrix(alpha, &h2->u->S, xt, yt);
    addevaltrans_amatrix(alpha, &h2->u->S, xta, yta);
  }
  else {
    assert(h2->son != 0);

    rsons = h2->rsons;
    csons = h2->csons;

    xtoff = cb->k;
    for (j = 0; j < csons; j++) {
      xt1 = 0;
      yta1 = 0;
      if (h2->son[j * rsons]->cb == cb) {
	xt1 = init_sub_avector(&tmp1, xt, cb->ktree, 0);
	yta1 = init_sub_avector(&tmp2, yta, cb->ktree, 0);
	xtoff += cb->t->nidx;
      }
      else {
	assert(j < cb->sons);

	xt1 = init_sub_avector(&tmp1, xt, cb->son[j]->ktree, xtoff);
	yta1 = init_sub_avector(&tmp2, yta, cb->son[j]->ktree, xtoff);
	xtoff += cb->son[j]->ktree;
      }

      ytoff = rb->k;
      for (i = 0; i < rsons; i++) {
	yt1 = 0;
	xta1 = 0;
	if (h2->son[i]->rb == rb) {
	  yt1 = init_sub_avector(&tmp3, yt, rb->ktree, 0);
	  xta1 = init_sub_avector(&tmp4, xta, rb->ktree, 0);
	  ytoff += rb->t->nidx;
	}
	else {
	  assert(i < rb->sons);

	  yt1 = init_sub_avector(&tmp3, yt, rb->son[i]->ktree, ytoff);
	  xta1 = init_sub_avector(&tmp4, xta, rb->son[i]->ktree, ytoff);
	  ytoff += rb->son[i]->ktree;
	}

	addevalsymm_offdiag(alpha, h2->son[i + j * rsons], xt1, xta1, yt1,
			    yta1);

	uninit_avector(xta1);
	uninit_avector(yt1);
      }
      assert(ytoff == rb->ktree);

      uninit_avector(yta1);
      uninit_avector(xt1);
    }
    assert(xtoff == cb->ktree);
  }
}

static void
addevalsymm_diag(field alpha, pch2matrix h2, pavector xt,
		 pavector xta, pavector yt, pavector yta)
{
  avector   tmp1, tmp2, tmp3, tmp4;
  pavector  xt1, xta1, yt1, yta1;
  pavector  xp, yp;
  pcclusterbasis rb = h2->rb;
  pcclusterbasis cb = h2->cb;
  pfield    aa;
  uint      lda, sons;
  uint      xtoff, ytoff;
  uint      n;
  uint      i, j;

  assert(h2->rb->t == h2->cb->t);
  assert(xt->size == h2->cb->ktree);
  assert(xta->size == h2->rb->ktree);
  assert(yt->size == h2->rb->ktree);
  assert(yta->size == h2->cb->ktree);

  if (h2->f) {
    aa = h2->f->a;
    lda = h2->f->ld;

    n = rb->t->nidx;
    xp = init_sub_avector(&tmp1, xt, n, cb->k);
    yp = init_sub_avector(&tmp2, yt, n, rb->k);

    for (j = 0; j < n; j++) {
      yp->v[j] += alpha * aa[j + j * lda] * xp->v[j];
      for (i = j + 1; i < n; i++) {
	yp->v[i] += alpha * aa[i + j * lda] * xp->v[j];
	yp->v[j] += alpha * CONJ(aa[i + j * lda]) * xp->v[i];
      }
    }

    uninit_avector(yp);
    uninit_avector(xp);
  }
  else {
    assert(h2->son != 0);
    assert(h2->rsons == h2->csons);

    sons = h2->rsons;

    xtoff = cb->k;
    for (j = 0; j < sons; j++) {
      xt1 = 0;
      yta1 = 0;
      if (h2->son[j * sons]->cb == cb) {
	xt1 = init_sub_avector(&tmp1, xt, cb->ktree, 0);
	yta1 = init_sub_avector(&tmp2, yta, cb->ktree, 0);
	xtoff += cb->t->nidx;
      }
      else {
	assert(j < cb->sons);

	xt1 = init_sub_avector(&tmp1, xt, cb->son[j]->ktree, xtoff);
	yta1 = init_sub_avector(&tmp2, yta, cb->son[j]->ktree, xtoff);
	xtoff += cb->son[j]->ktree;
      }

      ytoff = rb->k;

      yt1 = 0;
      xta1 = 0;
      if (h2->son[0]->rb == rb) {
	yt1 = init_sub_avector(&tmp3, yt, rb->ktree, 0);
	xta1 = init_sub_avector(&tmp4, xta, rb->ktree, 0);
	ytoff += rb->t->nidx;
      }
      else {
	/* Skip superdiagonal blocks */
	for (i = 0; i < j; i++)
	  ytoff += rb->son[i]->ktree;

	/* Subvectors for diagonal block */
	yt1 = init_sub_avector(&tmp3, yt, rb->son[i]->ktree, ytoff);
	xta1 = init_sub_avector(&tmp4, xta, rb->son[i]->ktree, ytoff);
	ytoff += rb->son[i]->ktree;
      }

      addevalsymm_diag(alpha, h2->son[j + j * sons], xt1, xta1, yt1, yta1);

      uninit_avector(xta1);
      uninit_avector(yt1);

      for (i = j + 1; i < sons; i++) {
	assert(i < rb->sons);

	yt1 = init_sub_avector(&tmp3, yt, rb->son[i]->ktree, ytoff);
	xta1 = init_sub_avector(&tmp4, xta, rb->son[i]->ktree, ytoff);
	ytoff += rb->son[i]->ktree;

	addevalsymm_offdiag(alpha, h2->son[i + j * sons], xt1, xta1, yt1,
			    yta1);

	uninit_avector(xta1);
	uninit_avector(yt1);
      }
      assert(ytoff == rb->ktree);

      uninit_avector(yta1);
      uninit_avector(xt1);
    }
    assert(xtoff == cb->ktree);
  }
}

void
addevalsymm_h2matrix(field alpha, pch2matrix h2, pcavector x, pavector y)
{
  pavector  xt, yt, xta, yta;

  assert(h2->rb->t == h2->cb->t);

  /* Transformed coefficients */
  xt = new_avector(h2->cb->ktree);
  xta = new_avector(h2->rb->ktree);
  yt = new_avector(h2->rb->ktree);
  yta = new_avector(h2->cb->ktree);

  /* Clear row coefficients */
  clear_avector(yt);
  clear_avector(yta);

  /* Column coefficients filled by forward transformation */
  forward_clusterbasis(h2->cb, x, xt);
  forward_clusterbasis(h2->rb, x, xta);

  /* Multiplication step */
  addevalsymm_diag(alpha, h2, xt, xta, yt, yta);

  /* Row coefficients added to result by backward transformation */
  backward_clusterbasis(h2->rb, yt, y);
  backward_clusterbasis(h2->cb, yta, y);

  /* Clean up */
  uninit_avector(yta);
  uninit_avector(yt);
  uninit_avector(xta);
  uninit_avector(xt);
}

/* ------------------------------------------------------------
 * On-the-fly matrix-vector multiplication
 * ------------------------------------------------------------ */

void
mvm_otf_h2matrix(field alpha, bool h2trans,
		 pcblock b, pcclusterbasis rb, pcclusterbasis cb,
		 buildN_t buildN, buildS_t buildS, void *data,
		 pcavector x, pavector y)
{
  if (h2trans)
    addevaltrans_otf_h2matrix(alpha, b, rb, cb, buildN, buildS, data, x, y);
  else
    addeval_otf_h2matrix(alpha, b, rb, cb, buildN, buildS, data, x, y);
}

void
fastaddeval_otf_h2matrix(field alpha,
			 pcblock b, pcclusterbasis rb, pcclusterbasis cb,
			 buildN_t buildN, buildS_t buildS, void *data,
			 pcavector xt, pavector yt)
{
  amatrix   loc1;
  pamatrix  N, S;
  avector   loc2, loc3;
  pavector  xp, yp, xt1, yt1;
  uint      rsons, csons;
  uint      xtoff, ytoff;
  uint      i, j;

  assert(xt->size == cb->ktree);
  assert(yt->size == rb->ktree);
  assert(b->rc == rb->t);
  assert(b->cc == cb->t);

  if (b->son) {
    rsons = b->rsons;
    csons = b->csons;

    if (b->son[0]->rc == b->rc) {
      assert(rsons == 1);

      if (b->son[0]->cc == b->cc) {
	assert(csons == 1);

	fastaddeval_otf_h2matrix(alpha, b->son[0], rb, cb,
				 buildN, buildS, data, xt, yt);
      }
      else {
	assert(csons == cb->sons);

	xtoff = cb->k;
	for (j = 0; j < csons; j++) {
	  xt1 =
	    init_sub_avector(&loc2, (pavector) xt, cb->son[j]->ktree, xtoff);

	  fastaddeval_otf_h2matrix(alpha, b->son[j], rb, cb->son[j],
				   buildN, buildS, data, xt1, yt);

	  uninit_avector(xt1);

	  xtoff += cb->son[j]->ktree;
	}
	assert(xtoff == cb->ktree);
      }
    }
    else {
      assert(rsons == rb->sons);

      if (b->son[0]->cc == b->cc) {
	assert(csons == 1);

	ytoff = rb->k;
	for (i = 0; i < rsons; i++) {
	  yt1 = init_sub_avector(&loc3, yt, rb->son[i]->ktree, ytoff);

	  fastaddeval_otf_h2matrix(alpha, b->son[i], rb->son[i], cb,
				   buildN, buildS, data, xt, yt1);

	  uninit_avector(yt1);

	  ytoff += rb->son[i]->ktree;
	}
	assert(ytoff == rb->ktree);
      }
      else {
	assert(csons == cb->sons);

	xtoff = cb->k;
	for (j = 0; j < csons; j++) {
	  xt1 =
	    init_sub_avector(&loc2, (pavector) xt, cb->son[j]->ktree, xtoff);

	  ytoff = rb->k;
	  for (i = 0; i < rsons; i++) {
	    yt1 = init_sub_avector(&loc3, yt, rb->son[i]->ktree, ytoff);

	    fastaddeval_otf_h2matrix(alpha, b->son[i + j * rsons], rb->son[i],
				     cb->son[j], buildN, buildS, data, xt1,
				     yt1);

	    uninit_avector(yt1);

	    ytoff += rb->son[i]->ktree;
	  }
	  assert(ytoff == rb->ktree);

	  uninit_avector(xt1);

	  xtoff += cb->son[j]->ktree;
	}
	assert(xtoff == cb->ktree);
      }
    }
  }
  else if (b->adm) {
    xt1 = init_sub_avector(&loc2, (pavector) xt, cb->k, 0);
    yt1 = init_sub_avector(&loc3, yt, rb->k, 0);

    S = init_amatrix(&loc1, rb->k, cb->k);
    buildS(data, rb->t, cb->t, S);

    addeval_amatrix(alpha, S, xt1, yt1);

    uninit_amatrix(S);
    uninit_avector(yt1);
    uninit_avector(xt1);
  }
  else {
    xp = init_sub_avector(&loc2, (pavector) xt, cb->t->nidx, cb->k);
    yp = init_sub_avector(&loc3, yt, rb->t->nidx, rb->k);

    N = init_amatrix(&loc1, rb->t->nidx, cb->t->nidx);
    buildN(data, rb->t->idx, cb->t->idx, N);

    addeval_amatrix(alpha, N, xp, yp);

    uninit_amatrix(N);
    uninit_avector(yp);
    uninit_avector(xp);
  }
}

void
addeval_otf_h2matrix(field alpha,
		     pcblock b, pcclusterbasis rb, pcclusterbasis cb,
		     buildN_t buildN, buildS_t buildS, void *data,
		     pcavector x, pavector y)
{
  pavector  xt, yt;

  xt = new_avector(cb->ktree);
  yt = new_avector(rb->ktree);

  clear_avector(yt);

  forward_clusterbasis(cb, x, xt);

  fastaddeval_otf_h2matrix(alpha, b, rb, cb, buildN, buildS, data, xt, yt);

  backward_clusterbasis(rb, yt, y);

  del_avector(yt);
  del_avector(xt);
}

void
fastaddevaltrans_otf_h2matrix(field alpha,
			      pcblock b, pcclusterbasis rb, pcclusterbasis cb,
			      buildN_t buildN, buildS_t buildS, void *data,
			      pcavector xt, pavector yt)
{
  amatrix   loc1;
  pamatrix  N, S;
  avector   loc2, loc3;
  pavector  xp, yp, xt1, yt1;
  uint      rsons, csons;
  uint      xtoff, ytoff;
  uint      i, j;

  assert(xt->size == rb->ktree);
  assert(yt->size == cb->ktree);
  assert(b->rc == rb->t);
  assert(b->cc == cb->t);

  if (b->son) {
    rsons = b->rsons;
    csons = b->csons;

    if (b->son[0]->rc == b->rc) {
      assert(rsons == 1);

      if (b->son[0]->cc == b->cc) {
	assert(csons == 1);

	fastaddevaltrans_otf_h2matrix(alpha, b->son[0], rb, cb,
				      buildN, buildS, data, xt, yt);
      }
      else {
	assert(csons == cb->sons);

	ytoff = cb->k;
	for (j = 0; j < csons; j++) {
	  yt1 = init_sub_avector(&loc2, yt, cb->son[j]->ktree, ytoff);

	  fastaddevaltrans_otf_h2matrix(alpha, b->son[j], rb, cb->son[j],
					buildN, buildS, data, xt, yt1);

	  uninit_avector(yt1);

	  ytoff += cb->son[j]->ktree;
	}
	assert(ytoff == cb->ktree);
      }
    }
    else {
      assert(rsons == rb->sons);

      if (b->son[0]->cc == b->cc) {
	assert(csons == 1);

	xtoff = rb->k;
	for (i = 0; i < rsons; i++) {
	  xt1 =
	    init_sub_avector(&loc3, (pavector) xt, rb->son[i]->ktree, xtoff);

	  fastaddevaltrans_otf_h2matrix(alpha, b->son[i], rb->son[i], cb,
					buildN, buildS, data, xt1, yt);

	  uninit_avector(xt1);

	  xtoff += rb->son[i]->ktree;
	}
	assert(xtoff == rb->ktree);
      }
      else {
	assert(csons == cb->sons);

	ytoff = cb->k;
	for (j = 0; j < csons; j++) {
	  yt1 = init_sub_avector(&loc2, yt, cb->son[j]->ktree, ytoff);

	  xtoff = rb->k;
	  for (i = 0; i < rsons; i++) {
	    xt1 =
	      init_sub_avector(&loc3, (pavector) xt, rb->son[i]->ktree,
			       xtoff);

	    fastaddevaltrans_otf_h2matrix(alpha, b->son[i + j * rsons],
					  rb->son[i], cb->son[j], buildN,
					  buildS, data, xt1, yt1);

	    uninit_avector(xt1);

	    xtoff += rb->son[i]->ktree;
	  }
	  assert(xtoff == rb->ktree);

	  uninit_avector(yt1);

	  ytoff += cb->son[j]->ktree;
	}
	assert(ytoff == cb->ktree);
      }
    }
  }
  else if (b->adm) {
    yt1 = init_sub_avector(&loc2, yt, cb->k, 0);
    xt1 = init_sub_avector(&loc3, (pavector) xt, rb->k, 0);

    S = init_amatrix(&loc1, rb->k, cb->k);
    buildS(data, rb->t, cb->t, S);

    addevaltrans_amatrix(alpha, S, xt1, yt1);

    uninit_amatrix(S);
    uninit_avector(yt1);
    uninit_avector(xt1);
  }
  else {
    yp = init_sub_avector(&loc2, yt, cb->t->nidx, cb->k);
    xp = init_sub_avector(&loc3, (pavector) xt, rb->t->nidx, rb->k);

    N = init_amatrix(&loc1, rb->t->nidx, cb->t->nidx);
    buildN(data, rb->t->idx, cb->t->idx, N);

    addevaltrans_amatrix(alpha, N, xp, yp);

    uninit_amatrix(N);
    uninit_avector(yp);
    uninit_avector(xp);
  }
}

void
addevaltrans_otf_h2matrix(field alpha,
			  pcblock b, pcclusterbasis rb, pcclusterbasis cb,
			  buildN_t buildN, buildS_t buildS, void *data,
			  pcavector x, pavector y)
{
  pavector  xt, yt;

  xt = new_avector(rb->ktree);
  yt = new_avector(cb->ktree);

  clear_avector(yt);

  forward_clusterbasis(rb, x, xt);

  fastaddevaltrans_otf_h2matrix(alpha, b, rb, cb, buildN, buildS, data, xt,
				yt);

  backward_clusterbasis(cb, yt, y);

  del_avector(yt);
  del_avector(xt);
}

/* ------------------------------------------------------------
 * Spectral norm
 * ------------------------------------------------------------ */

real
norm2_h2matrix(pch2matrix H2)
{
  return norm2_matrix((mvm_t) mvm_h2matrix, (void *) H2,
		      H2->rb->t->nidx, H2->cb->t->nidx);
}

real
norm2diff_h2matrix(pch2matrix a, pch2matrix b)
{
  return norm2diff_matrix((mvm_t) mvm_h2matrix, (void *) a,
			  (mvm_t) mvm_h2matrix, (void *) b,
			  a->rb->t->nidx, a->cb->t->nidx);
}

real
norm2diff_amatrix_h2matrix(pch2matrix a, pcamatrix b)
{
  return norm2diff_matrix((mvm_t) mvm_h2matrix, (void *) a,
			  (mvm_t) mvm_amatrix, (void *) b,
			  a->rb->t->nidx, a->cb->t->nidx);
}

/* ------------------------------------------------------------
 * Drawing
 * ------------------------------------------------------------ */

#ifdef USE_CAIRO
static void
cairodraw(cairo_t * cr, pch2matrix G, bool storage, uint levels)
{
  cairo_text_extents_t extents;
  char      buf[11];
  uint      rsons, csons;
  uint      rsize, csize;
  uint      roff, coff;
  uint      i, j;

  if (G->son && levels != 1) {
    rsons = G->rsons;
    csons = G->csons;

    coff = 0;
    for (j = 0; j < csons; j++) {
      roff = 0;
      for (i = 0; i < rsons; i++) {
	cairo_save(cr);
	cairo_translate(cr, coff, roff);
	cairodraw(cr, G->son[i + j * rsons], storage, levels - 1);
	cairo_restore(cr);

	roff += G->son[i + j * rsons]->rb->t->nidx;
      }
      assert(roff == G->rb->t->nidx);

      coff += G->son[j * rsons]->cb->t->nidx;
    }
    assert(coff == G->cb->t->nidx);
  }
  else {
    rsize = G->rb->t->nidx;
    csize = G->cb->t->nidx;

    if (G->son) {
      cairo_rectangle(cr, 0.0, 0.0, csize, rsize);
      cairo_save(cr);
      cairo_set_source_rgb(cr, 0.9, 0.9, 1.0);
      cairo_fill_preserve(cr);
      cairo_restore(cr);
      cairo_stroke(cr);
    }
    else if (G->u) {
      if (storage) {
	cairo_rectangle(cr, 0.0, 0.0, (G->cb->k > csize ? csize : G->cb->k),
			(G->rb->k > rsize ? rsize : G->rb->k));
	cairo_save(cr);
	cairo_set_source_rgb(cr, 1.0, 0.0, 1.0);
	cairo_fill(cr);
	cairo_restore(cr);

	if (G->rb->k > 0 && G->cb->k > 0) {
	  sprintf(buf, "%dx%d", G->rb->k, G->cb->k);
	  cairo_set_font_size(cr, UINT_MIN(rsize, csize) / 4.75);
	  cairo_text_extents(cr, buf, &extents);

	  cairo_save(cr);
	  cairo_translate(cr, (csize - extents.width) / 2,
			  (rsize + extents.height) / 2);
	  cairo_show_text(cr, buf);
	  cairo_restore(cr);
	}

	cairo_rectangle(cr, 0.0, 0.0, csize, rsize);
	cairo_stroke(cr);
      }
      else {
	cairo_rectangle(cr, 0.0, 0.0, csize, rsize);
	cairo_save(cr);
	cairo_set_source_rgb(cr, 0.2, 0.2, 1.0);
	cairo_fill_preserve(cr);
	cairo_restore(cr);
	cairo_stroke(cr);
      }
    }
    else if (G->f) {
      cairo_rectangle(cr, 0.0, 0.0, csize, rsize);
      cairo_save(cr);
      cairo_set_source_rgb(cr, 1.0, 0.0, 0.0);
      cairo_fill_preserve(cr);
      cairo_restore(cr);
      cairo_stroke(cr);
    }
    else {
      cairo_rectangle(cr, 0.0, 0.0, csize, rsize);
      cairo_stroke(cr);
    }
  }
}

void
draw_cairo_h2matrix(cairo_t * cr, pch2matrix G, bool storage, uint levels)
{
  double    sx, sy, ex, ey;
  uint      rsize, csize;
  double    scalex, scaley, scale;

  /* Save Cairo state */
  cairo_save(cr);

  /* Obtain size of block */
  rsize = G->rb->t->nidx;
  csize = G->cb->t->nidx;

  /* Obtain size of current Cairo bounding box */
  cairo_clip_extents(cr, &sx, &sy, &ex, &ey);

  /* Compute scaling factor */
  scalex = (ex - sx) / rsize;
  scaley = (ey - sy) / csize;
  scale = (scalex < scaley ? scalex : scaley);

  /* Center block in bounding box */
  cairo_translate(cr, 0.5 * (ex - sx - scale * rsize),
		  0.5 * (ey - sy - scale * csize));

  /* Scale coordinates */
  cairo_scale(cr, scale, scale);
  cairo_set_line_width(cr, cairo_get_line_width(cr) / scale);

  /* Start drawing */
  cairodraw(cr, G, storage, levels);

  /* Restore Cairo state */
  cairo_restore(cr);
}
#endif
