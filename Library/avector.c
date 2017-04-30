
/* ------------------------------------------------------------
 * This is the file "avector.c" of the KIPS package.
 * All rights reserved, Steffen Boerm 2009
 * ------------------------------------------------------------ */

#include "avector.h"

#include "settings.h"
#include "basic.h"

#include <math.h>
#include <stdio.h>

static uint active_avector = 0;

/* ------------------------------------------------------------
 * Constructors and destructors
 * ------------------------------------------------------------ */

pavector
init_avector(pavector v, uint size)
{
  assert(v != NULL);

  v->v = (size > 0 ? allocfield(size) : NULL);
  v->size = size;
  v->owner = NULL;

#ifdef USE_OPENMP
#pragma omp atomic
#endif
  active_avector++;

  return v;
}

pavector
init_sub_avector(pavector v, pavector src, uint size, uint off)
{
  assert(v != NULL);
  assert(src != NULL);
  assert(off + size <= src->size);

  v->v = src->v + off;
  v->size = size;
  v->owner = src;

#ifdef USE_OPENMP
#pragma omp atomic
#endif
  active_avector++;

  return v;
}

pavector
init_column_avector(pavector v, pamatrix src, uint col)
{
  longindex lda;

  assert(v != NULL);
  assert(src != NULL);
  assert(col < src->cols);

  lda = src->ld;

  v->v = src->a + col * lda;
  v->size = src->rows;
  v->owner = src;

#ifdef USE_OPENMP
#pragma omp atomic
#endif
  active_avector++;

  return v;
}

void
uninit_avector(pavector v)
{
  if (!v->owner)
    freemem(v->v);

  assert(active_avector > 0);

#ifdef USE_OPENMP
#pragma omp atomic
#endif
  active_avector--;
}

pavector
new_avector(uint size)
{
  pavector  v;

  v = (pavector) allocmem(sizeof(avector));

  init_avector(v, size);

  return v;
}

pavector
new_sub_avector(pavector src, uint size, uint off)
{
  pavector  v;

  v = (pavector) allocmem(sizeof(avector));

  init_sub_avector(v, src, size, off);

  return v;
}

void
del_avector(pavector v)
{
  uninit_avector(v);
  freemem(v);
}

void
resize_avector(pavector v, uint size)
{
  assert(v->owner == NULL);

  if (size != v->size) {
    freemem(v->v);
    v->v = allocfield(size);
    v->size = size;
  }
}

/* ------------------------------------------------------------
 * Statistics
 * ------------------------------------------------------------ */

uint
getactives_avector()
{
  return active_avector;
}

size_t
getsize_avector(pcavector v)
{
  size_t    sz;

  sz = sizeof(avector);
  if (v->owner == NULL)
    sz += (size_t) sizeof(field) * v->size;

  return sz;
}

size_t
getsize_heap_avector(pcavector v)
{
  size_t    sz;

  sz = 0;
  if (v->owner == NULL)
    sz += (size_t) sizeof(field) * v->size;

  return sz;
}

/* ------------------------------------------------------------
 * Simple utility functions
 * ------------------------------------------------------------ */

void
clear_avector(pavector v)
{
  uint      i;

  for (i = 0; i < v->size; i++)
    v->v[i] = 0.0;
}

void
fill_avector(field alpha, pavector v)
{
  uint      i;

  for (i = 0; i < v->size; i++)
    v->v[i] = alpha;
}

void
random_avector(pavector v)
{
  uint      i;

  for (i = 0; i < v->size; i++) {
    v->v[i] = FIELD_RAND();
  }
}

void
copy_avector(pcavector v, pavector w)
{
  uint      i;

  assert(v->size == w->size);

  for (i = 0; i < v->size; i++)
    w->v[i] = v->v[i];
}

void
copy_sub_avector(pcavector v, pavector w)
{
  uint      i, n;

  n = UINT_MIN(v->size, w->size);

  for (i = 0; i < n; i++)
    w->v[i] = v->v[i];
}

pavector
clone_avector(pcavector src)
{
  pavector x;

  x = new_avector(src->size);
  copy_avector(src, x);

  return x;
}

void
print_avector(pcavector v)
{
  uint      size = v->size;
  uint      i;

  (void) printf("avector(%u)\n", size);
  if (size == 0)
    return;

#ifdef USE_COMPLEX
  (void) printf("  (%f+%fi", REAL(v->v[0]), IMAG(v->v[0]));
  for (i = 1; i < size; i++)
    (void) printf(" %f+%fi", REAL(v->v[i]), IMAG(v->v[i]));
  (void) printf(")\n");
#else
  (void) printf("  (%f", v->v[0]);
  for (i = 1; i < size; i++)
    (void) printf(" %f", v->v[i]);
  (void) printf(")\n");
#endif
}

/* ------------------------------------------------------------
 * Very basic linear algebra
 * ------------------------------------------------------------ */

void
scale_avector(field alpha, pavector v)
{
  scal(v->size, alpha, v->v, 1);
}

real
norm2_avector(pcavector v)
{
  return nrm2(v->size, v->v, 1);
}

field
dotprod_avector(pcavector x, pcavector y)
{
  assert(x->size == y->size);

  return dot(true, false, x->size, x->v, 1, y->v, 1);
}

void
add_avector(field alpha, pcavector x, pavector y)
{
  assert(y->size == x->size);

  axpy(x->size, alpha, x->v, 1, y->v, 1);
}
