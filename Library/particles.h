
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

struct _particles {
  /** @brief Number of particles */
  uint n;

  /** @brief First coordinates */
  real *x0;

  /** @brief Second coordinates */
  real *x1;

  /** @brief Third coordinates */
  real *x2;
};

HEADER_PREFIX pparticles
new_particles(uint n);

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

/** @} */

#endif
