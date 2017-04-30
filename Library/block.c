
/* ------------------------------------------------------------
 * This is the file "block.c" of the KIPS package.
 * All rights reserved, Steffen Boerm 2010
 * ------------------------------------------------------------ */

#include "block.h"

#include "basic.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>

static uint active_block = 0;

/* ------------------------------------------------------------
 * Constructors and destructors
 * ------------------------------------------------------------ */

pblock
new_block(pcluster rc, pcluster cc, uint rsons, uint csons)
{
  pblock b;
  uint i, j;

  b = (pblock) allocmem(sizeof(block));
  b->rc = rc;
  b->cc = cc;
  b->adm = false;
  b->rsons = rsons;
  b->csons = csons;
  b->desc = 0;

  b->son = NULL;
  if(rsons > 0 && csons > 0) {
    b->son = (pblock *) allocmem((size_t) sizeof(block) * rsons * csons);
    for(j=0; j<csons; j++)
      for(i=0; i<rsons; i++)
	b->son[i+j*rsons] = NULL;
  }

#ifdef USE_OPENMP
#pragma omp atomic
#endif
  active_block++;

  return b;
}

void
update_block(pblock b)
{
  uint desc;
  uint i, j;

  desc = 1;
  for(j=0; j<b->csons; j++)
    for(i=0; i<b->rsons; i++)
      desc += b->son[i+j*b->rsons]->desc;

  b->desc = desc;
}

void
del_block(pblock b)
{
  uint i, j;

  if(b->son) {
    for(j=0; j<b->csons; j++)
      for(i=0; i<b->rsons; i++)
	del_block(b->son[i+j*b->rsons]);
    freemem(b->son);
  }

  freemem(b);

  assert(active_block > 0);

#ifdef USE_OPENMP
#pragma omp atomic
#endif
  active_block--;
}

/* ------------------------------------------------------------
 * Statistics
 * ------------------------------------------------------------ */

uint
getactives_block()
{
  return active_block;
}

size_t
getsize_block(pcblock b)
{
  size_t sz;
  uint i, j;

  sz = (size_t) sizeof(block);
  if(b->son) {
    sz += (size_t) sizeof(pblock) * b->rsons * b->csons;
    for(j=0; j<b->csons; j++)
      for(i=0; i<b->rsons; i++)
	sz += getsize_block(b->son[i+j*b->rsons]);
  }

  return sz;
}

uint
getdepth_block(pcblock b)
{
  uint depth, depth1;
  uint i, j;

  depth = 0;
  if(b->son) {
    for(j=0; j<b->csons; j++)
      for(i=0; i<b->rsons; i++) {
	depth1 = getdepth_block(b->son[i+j*b->rsons]) + 1;
	if(depth1 > depth)
	  depth = depth1;
      }
  }
  return depth;
}

/* ------------------------------------------------------------
 * Drawing
 * ------------------------------------------------------------ */

#ifdef USE_CAIRO
static void
cairodraw_subblock(cairo_t *cr, pcblock b, int levels)
{
  uint rsons, csons;
  uint rsize, csize;
  uint roff, coff;
  uint i, j;

  if(b->son && levels != 1) {
    rsons = b->rsons;
    csons = b->csons;

    coff = 0;
    for(j=0; j<csons; j++) {
      roff = 0;
      for(i=0; i<rsons; i++) {
	cairo_save(cr);
	cairo_translate(cr, coff, roff);
	cairodraw_subblock(cr, b->son[i+j*rsons], levels-1);
	cairo_restore(cr);

	roff += b->son[i+j*rsons]->rc->size;
      }
      assert(roff == b->rc->size);

      coff += b->son[j*rsons]->cc->size;
    }
    assert(coff == b->cc->size);
  }
  else {
    rsize = b->rc->size;
    csize = b->cc->size;
    
    if(b->son) {
      cairo_rectangle(cr, 0.0, 0.0, csize, rsize);
      cairo_stroke(cr);
    }
    else if(b->adm) {
      cairo_rectangle(cr, 0.0, 0.0, csize, rsize);
      cairo_save(cr);
      cairo_set_source_rgb(cr, 0.2, 0.2, 1.0);
      cairo_fill_preserve(cr);
      cairo_restore(cr);
      cairo_stroke(cr);
    }
    else {
      cairo_rectangle(cr, 0.0, 0.0, csize, rsize);
      cairo_save(cr);
      cairo_set_source_rgb(cr, 1.0, 0.0, 0.0);
      cairo_fill_preserve(cr);
      cairo_restore(cr);
      cairo_stroke(cr);
    }
  }
}

void
cairodraw_block(cairo_t *cr, pcblock b, int levels)
{
  double sx, sy, ex, ey;
  uint rsize, csize;
  double scalex, scaley, scale;

  /* Save Cairo state */
  cairo_save(cr);

  /* Obtain size of block */
  rsize = b->rc->size;
  csize = b->cc->size;

  /* Obtain size of current Cairo bounding box */
  cairo_clip_extents(cr, &sx, &sy, &ex, &ey);

  /* Compute scaling factor */
  scalex = (ex - sx) / rsize;
  scaley = (ey - sy) / csize;
  scale = (scalex < scaley ? scalex : scaley);

  /* Center block in bounding box */
  cairo_translate(cr,
		  0.5 * (ex - sx - scale * rsize),
		  0.5 * (ey - sy - scale * csize));

  /* Scale coordinates */
  cairo_scale(cr, scale, scale);
  cairo_set_line_width(cr, cairo_get_line_width(cr) / scale);

  /* Start drawing */
  cairodraw_subblock(cr, b, levels);

  /* Restore Cairo state */
  cairo_restore(cr);
}
#endif

/* ------------------------------------------------------------
 * Building standard block trees
 * ------------------------------------------------------------ */

pblock
build_block(pcluster rc, pcluster cc,
	    bool (*admissible)(pcluster rc, pcluster cc, void *data),
	    void *data,
	    uint levels, bool strict)
{
  pblock b;
  uint i, j;

  b = 0;
  if(admissible(rc, cc, data)) {
    b = new_block(rc, cc, 0, 0);
    b->adm = true;
  }
  else if((rc->sons == 0 && cc->sons == 0) ||
	  (!strict && (rc->sons == 0 || cc->sons == 0)) ||
	  levels <= 0) {
    b = new_block(rc, cc, 0, 0);
  }
  else {
    if(rc->sons > 0) {
      if(cc->sons > 0) {
	b = new_block(rc, cc, rc->sons, cc->sons);
	for(j=0; j<cc->sons; j++)
	  for(i=0; i<rc->sons; i++)
	    b->son[i+j*rc->sons] = build_block(rc->son[i], cc->son[j],
					       admissible, data,
					       levels-1, strict);
      }
      else {
	b = new_block(rc, cc, rc->sons, 1);
	for(i=0; i<rc->sons; i++)
	  b->son[i] = build_block(rc->son[i], cc,
				  admissible, data,
				  levels-1, strict);
      }
    }
    else {
      if(cc->sons > 0) {
	b = new_block(rc, cc, 1, cc->sons);
	for(j=0; j<cc->sons; j++)
	  b->son[j] = build_block(rc, cc->son[j],
				  admissible, data,
				  levels-1, strict);
      }
      else {
	/* Should be impossible */
	assert(0);
      }
    }
  }

  update_block(b);

  return b;
}

bool
h2std_admissibility(pcluster rc, pcluster cc, void *data)
{
  preal pe = (preal) data;
  real eta = (*pe);
  real rdiam, cdiam, dist;

  rdiam = diam_cluster(rc);
  cdiam = diam_cluster(cc);
  dist = 2.0 * eta * dist_cluster(rc, cc);

  return (rdiam <= dist && cdiam <= dist);
}

pblock
buildh2std_block(pcluster rc, pcluster cc, real eta)
{
  return build_block(rc, cc, h2std_admissibility, &eta, 64, true);
}

bool
hstd_admissibility(pcluster rc, pcluster cc, void *data)
{
  preal pe = (preal) data;
  real rdiam, cdiam, dist;

  rdiam = diam_cluster(rc);
  cdiam = diam_cluster(cc);
  dist = 2.0 * (*pe) * dist_cluster(rc, cc);

  return (rdiam <= dist || cdiam <= dist);
}

pblock
buildhstd_block(pcluster rc, pcluster cc, real eta)
{
  return build_block(rc, cc, hstd_admissibility, &eta, 64, false);
}
