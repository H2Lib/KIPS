
/* ------------------------------------------------------------
 * This is the file "particles.h" of the KIPS package.
 * All rights reserved, Steffen Boerm 2017
 * ------------------------------------------------------------ */

/** @file particles.h
 *  @author Steffen B&ouml;rm */

#ifndef PARTICLES_H
#define PARTICLES_H

/** @defgroup particles particles
 *  @brief Collection of particles.
 *  @{ */

typedef struct _particles particles;
typedef particles *pparticles;
typedef const particles *pcparticles;

#include "settings.h"
#include "avector.h"
#include "clustergeometry.h"
#include "interpolation.h"

struct _particles {
  /** @brief Number of particles */
  uint n;

  /** @brief Particle coordinates
   *
   *  <tt>x[i][j]</tt> is the <tt>j</tt>-th coordinate of the
   *  <tt>i</tt>-th particle */
  real **x;

  /** @brief <tt>x</tt> in row-major order */
  real *xdata;

  /** @brief Interpolation points */
  pinterpolation in;
};

HEADER_PREFIX pparticles
new_particles(uint n, uint m);

HEADER_PREFIX void
del_particles(pparticles p);

HEADER_PREFIX void
random_particles(pparticles p);

HEADER_PREFIX pclustergeometry
buildgeometry_particles(pcparticles p);

HEADER_PREFIX field
potential_newton(const real *x, const real *y);

HEADER_PREFIX void
addeval_direct_particles(field alpha, pcparticles p,
			 pcavector m, pavector phi);

HEADER_PREFIX void
buildV_particles(pcparticles p, pccluster t, pamatrix V);

HEADER_PREFIX void
buildE_particles(pcparticles p, pccluster son, pccluster father,
		 pamatrix E);

/** @} */

#endif
