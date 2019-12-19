/** ------------------------------------------------------------
 *    This is the file "lj.h" of the KIPS package.
 *    All rights reserved, Jonas Lorenzen 2019
 *  ------------------------------------------------------------ */

/** @file lj.h
 *  @author Jonas Lorenzen
 */

#ifndef LJ_H
#define LJ_H

/** @defgroup lj lj
 *  @brief Evaluation of the Lennard-Jones potential (and forces) in a periodic grid.
 *
 *  @attention All functions in this module enforce the use of SI units!
 *  @{ */

#include "basic.h"
#include "spatialgeometry.h"
#include "kernelmatrix.h"

/** @brief Representation of the necessary parameters for a kernel function 
 *  that evaluates a Lennard-Jones potential or forces. */
struct _lj {
  /** @brief Depth of the potential well. */
  real eps;
  
  /** @brief Zero of the potential to the power of 6. */
  real sig6;
  
  /** @brief Cutoff distance. */
  real R;
  
  /** @brief Information about the unit cell. */
  pcspatialgeometry sg;
  
  /** @brief Derived quantities for efficient evaluation of the potential. */
  real a;
  real b;
  real sig12;
};

/** @brief Type definition of Lennard-Jones parameters. */
typedef struct _lj lj;

/** @brief Pointer to a @ref lj object. */
typedef lj *plj;

/** ----------------------------------------------------------------
 *    Parameters
 *  ---------------------------------------------------------------- */

/** @brief Collection of necessary parameters for evaluation of a Lennard-Jones-type 
 *  kernel or gradient function by a smooth cutoff scheme.
 *
 *  Incorporates a shift of the Lennard-Jones potential to guarantee a smooth cutoff.
 *
 *  @param R Cutoff radius.
 *  @param sig Zero of the Lennard-Jones potential.
 *  @param eps Depth of the potential well.
 *  @param sg Information about the unit cell.
 *
 *  @returns Pointer to new @ref lj object with the given parameters. */
HEADER_PREFIX plj
setParameters_lj (real R, real sig, real eps, pcspatialgeometry sg);

/** @brief Collection of necessary parameters for evaluation of the short-range part 
 *  of a Lennard-Jones-type kernel or gradient function.
 *
 *  Uses a standard cutoff, which will lead to a discontinuous potential.
 *
 *  @param R Cutoff radius.
 *  @param sig Zero of the Lennard-Jones potential.
 *  @param eps Depth of the potential well.
 *  @param sg Information about the unit cell.
 *
 *  @returns Pointer to new @ref lj object with the given parameters. */
HEADER_PREFIX plj
setParametersUnshifted_lj (real R, real sig, real eps, pcspatialgeometry sg);

/** @brief Delete a set of Lennard-Jones parameters.
 *
 *  @param The given @ref lj object to be deleted. */
HEADER_PREFIX void 
delParameters_lj (plj lj);


/** ----------------------------------------------------------------
 *    Lennard-Jones potential
 *  ---------------------------------------------------------------- */

/** @brief Kernel function for the short-range part of the Lennard-Jones pair potential.
 *
 *  @param x Evaluation point.
 *  @param y Evaluation point.
 *  @param xmol Molecule number of @p x.
 *  @param ymol Molecule number of @p y.
 *  @param data Parameters of the Coulomb potential.
 *
 *  @returns Kernel function for Lennard-Jones pair potential, 
 *          evaluated in @p x and @p y. */
HEADER_PREFIX real
kernel_lj (pcreal x, pcreal y, uint xmol, uint ymol, void *data);

/** @brief Potential energy of a given system of Lennard-Jones sites.
 *
 *  @remark The input potential matrix is supposed to have the correct 
 *          structure, but leaf and nearfield matrices will be updated 
 *          within this function.
 * 
 *  @param klj Kernel matrix of the Lennard-Jones potential.
 *  @param sg Geometric information.
 *  @param Vlj Potential matrix.
 *
 *  @returns Value of the potential energy. */
HEADER_PREFIX real
energy_lj (pckernelmatrix klj, pspatialgeometry sg, ph2matrix Vlj);

/** ----------------------------------------------------------------
 *    Lennard-Jones forces
 *  ---------------------------------------------------------------- */

/** @brief Gradient function for the short-range part of the Lennard-Jones pair potential.
 *
 *  @param x Point charge location.
 *  @param y Point that exerts force.
 *  @param xmol Molecule number of @p x.
 *  @param ymol Molecule number of @p y.
 *  @param data Parameters of the Coulomb potential.
 *  @param f Gradient of the Coulombic pair potential, to be set by this function. */
HEADER_PREFIX void
gradient_lj (pcreal x, pcreal y, uint xmol, uint ymol, void *data, pfield f);

/** @brief Compute a concatenated vector of the forces acting on the Lennard-Jones 
 *        sites of the given system.
 *
 *  The i-th to (i+dim-1)-th entries of the vector correspond to the force vector 
 *  acting on the i-th site.
 *
 *  @remark The input gradient matrix is supposed to have the correct 
 *          structure, but leaf and nearfield matrices will be updated 
 *          within this function.
 *
 *  @param klj Kernel matrix for the Lennard-Jones potential.
 *  @param sg Geometric information.
 *  @param Flj Gradient matrix. 
 *  @param f Pointer to force vector. */
HEADER_PREFIX void
fullForce_lj (pckernelmatrix klj, pspatialgeometry sg, ph2matrix Flj, pavector f);

/** @} */

#endif
