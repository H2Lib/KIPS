
/* ------------------------------------------------------------
 * This is the file "interpolation.h" of the KIPS package.
 * All rights reserved, Steffen Boerm 2017
 * ------------------------------------------------------------ */

/** @file interpolation.h
 *  @author Steffen B&ouml;rm */

#ifndef INTERPOLATION_H
#define INTERPOLATION_H

/** @defgroup interpolation interpolation
 *  @brief One-dimensional Lagrange interpolation.
 *  @{ */

typedef struct _interpolation interpolation;
typedef interpolation *pinterpolation;
typedef const interpolation *pcinterpolation;

#include "settings.h"

struct _interpolation {
  /** @brief Number of interpolation points */
  uint m;

  /** @brief Left point of interval */
  real a;

  /** @brief Right point of interval */
  real b;

  /** @brief Interpolation points */
  real *xi;
};

HEADER_PREFIX pinterpolation
new_interpolation(uint m);

HEADER_PREFIX void
del_interpolation(pinterpolation in);

HEADER_PREFIX void
chebyshev_interpolation(real a, real b, pinterpolation in);

HEADER_PREFIX pinterpolation
build_transformed_interpolation(pcinterpolation in, real a, real b);

#ifdef __GNUC__
INLINE_PREFIX real
evallagrange_interpolation(pcinterpolation, uint, real) __attribute__((unused,const));
#endif

INLINE_PREFIX real
evallagrange_interpolation(pcinterpolation in, uint nu, real x)
{
  const real *xi = in->xi;
  uint m = in->m;
  real num, den;
  uint i;

  num = den = 1.0;
  for(i=0; i<nu; i++) {
    num *= x - xi[i];
    den *= xi[nu] - xi[i];
  }
  for(i=nu+1; i<m; i++) {
    num *= x - xi[i];
    den *= xi[nu] - xi[i];
  }

  return num / den;
}

#endif
