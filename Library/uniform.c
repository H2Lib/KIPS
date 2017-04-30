
/* ------------------------------------------------------------
 * This is the file "uniform.c" of the KIPS package.
 * All rights reserved, Steffen Boerm 2009
 * ------------------------------------------------------------ */

#include <math.h>

#include "uniform.h"
#include "basic.h"

/* ------------------------------------------------------------
 * Constructors and destructors
 * ------------------------------------------------------------ */

puniform
new_uniform(pclusterbasis rb, pclusterbasis cb)
{
  puniform  u;

  u = allocmem(sizeof(uniform));

  u->rb = rb;
  u->cb = cb;

  init_amatrix(&u->S, rb->k, cb->k);

  return u;
}

void
del_uniform(puniform u)
{
  assert(u != 0);

  uninit_amatrix(&u->S);

  u->rb = 0;
  u->cb = 0;

  freemem(u);
}

/* ------------------------------------------------------------
 * Statistics
 * ------------------------------------------------------------ */

size_t
getsize_uniform(pcuniform u)
{
  size_t    sz;

  sz = sizeof(uniform);
  sz += getsize_heap_amatrix(&u->S);

  return sz;
}

/* ------------------------------------------------------------
 * Simple utility functions
 * ------------------------------------------------------------ */

void
clear_uniform(puniform u)
{
  clear_amatrix(&u->S);
}

void
copy_uniform(bool trans, pcuniform src, puniform trg)
{
  if (trans) {
    assert(trg->rb == src->cb);
    assert(trg->cb == src->rb);

    copy_amatrix(true, &src->S, &trg->S);
  }
  else {
    assert(trg->rb == src->rb);
    assert(trg->cb == src->cb);

    copy_amatrix(false, &src->S, &trg->S);
  }
}

puniform
clone_uniform(pcuniform src)
{
  puniform  u;

  u = new_uniform(src->rb, src->cb);
  copy_uniform(false, src, u);

  return u;
}

void
scale_uniform(field alpha, puniform u)
{
  scale_amatrix(alpha, &u->S);
}

void
random_uniform(puniform u)
{
  random_amatrix(&u->S);
}
