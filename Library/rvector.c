
/* ------------------------------------------------------------
 * This is the file "rvector.c" of the KIPS package.
 * All rights reserved, Steffen Boerm 2015
 * ------------------------------------------------------------ */

#include "rvector.h"

#include "basic.h"

#include <math.h>
#include <stdio.h>

static uint active_rvector = 0;

/* ------------------------------------------------------------
   Constructors and destructors
   ------------------------------------------------------------ */

prvector
init_rvector(prvector v, uint size)
{
  assert(v != NULL);

  v->v = (size > 0 ? allocreal(size) : NULL);
  v->size = size;

#ifdef USE_OPENMP
#pragma omp atomic
#endif
  active_rvector++;

  return v;
}

void
uninit_rvector(prvector v)
{
  freemem(v->v);

  assert(active_rvector > 0);

#ifdef USE_OPENMP
#pragma omp atomic
#endif
  active_rvector--;
}

prvector
new_rvector(uint size)
{
  prvector v;

  v = (prvector) allocmem(sizeof(rvector));

  init_rvector(v, size);

  return v;
}

void
del_rvector(prvector v)
{
  uninit_rvector(v);
  freemem(v);
}

void
resize_rvector(prvector v, uint size)
{
  if(size != v->size) {
    freemem(v->v);
    v->v = allocreal(size);
    v->size = size;
  }
}

/* ------------------------------------------------------------
   Statistics
   ------------------------------------------------------------ */

uint
getactives_rvector()
{
  return active_rvector;
}

size_t
getsize_rvector(pcrvector v)
{
  size_t sz;

  sz = sizeof(rvector);
  sz += (size_t) sizeof(real) * v->size;

  return sz;
}

/* ------------------------------------------------------------
   Simple utility functions
   ------------------------------------------------------------ */

void
copy_rvector(pcrvector v, prvector w)
{
  uint i;

  assert(v->size == w->size);

  for(i=0; i<v->size; i++)
    w->v[i] = v->v[i];
}

void
print_rvector(pcrvector v)
{
  uint size = v->size;
  uint i;

  (void) printf("rvector(%u)\n", size);
  if(size == 0)
    return;

  (void) printf("  (% .5e", v->v[0]);
  for(i=1; i<size; i++)
    (void) printf(" % .5e", v->v[i]);
  (void) printf(")\n");
}
