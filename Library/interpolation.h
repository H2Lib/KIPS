
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

#endif
