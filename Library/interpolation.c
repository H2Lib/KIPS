
/* ------------------------------------------------------------
 * This is the file "interpolation.c" of the KIPS package.
 * All rights reserved, Steffen Boerm 2017
 * ------------------------------------------------------------ */

#include "interpolation.h"
#include "basic.h"

pinterpolation
new_interpolation(uint m)
{
  pinterpolation in;

  in = (pinterpolation) allocmem(sizeof(interpolation));
  in->m = m;
  in->a = 0.0;
  in->b = 0.0;
  in->xi = allocreal(m);

  return in;
}

void
del_interpolation(pinterpolation in)
{
  freemem(in->xi);
  in->xi = 0;
  freemem(in);
}

void
chebyshev_interpolation(real a, real b, pinterpolation in)
{
  preal xi = in->xi;
  uint m = in->m;
  real rad, mid;
  uint i;

  assert(a < b);

  in->a = a;
  in->b = b;

  rad = 0.5 * (b - a);
  mid = 0.5 * (b + a);

  for(i=0; i<m; i++)
    xi[i] = mid + rad * REAL_COS(M_PI * (i+0.5) / m);
}

pinterpolation
build_transformed_interpolation(pcinterpolation in, real a, real b)
{
  pinterpolation in2;
  pcreal xi = in->xi;
  preal xi2;
  uint m = in->m;
  real c0, c1;
  uint i;

  assert(a < b);
  assert(in->a < in->b);

  in2 = new_interpolation(m);
  in2->a = a;
  in2->b = b;
  xi2 = in2->xi;

  c0 = (in->b * a - in->a * b) / (in->b - in->a);
  c1 = (b - a) / (in->b - in->a);

  for(i=0; i<m; i++)
    xi2[i] = c0 + c1 * xi[i];

  return in2;
}
