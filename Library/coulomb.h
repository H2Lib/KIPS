/* ------------------------------------------------------------
 * This is the file "coulomb.h" of the KIPS package.
 * All rights reserved, Jonas Lorenzen 2018
 * ------------------------------------------------------------ */

/** @file coulomb.h
 *  @author Jonas Lorenzen
 */

#ifndef COULOMB_H
#define COULOMB_H

/** @defgroup coulomb coulomb
 *  @brief Evaluation of the Coulomb potential (and forces) in a periodic grid.
 *  @{ */

#include "basic.h"
#include "spatialgeometry.h"
#include "kernelmatrix.h"

/** @brief Splitting function for short and long range parts of the coulomb potential.
  This is a weight function for point distances up to a certain Cutoff radius.*/
typedef real (*split) (real alpha, real R, real r);

/** @brief Representation of the necessary parameters for a kernel function
  that evaluates a coulombic potential or forces. */
struct _parameters {
  /** @brief Range splitting function. */
  split s;
  
  /** @brief Derivative of the splitting function. */
  split sd;
  
  /** @brief Damping parameter for certain splitting functions. */
  real alpha;
  
  /** @brief Cutoff distance. */
  real R;
  
  /** @brief Information about the unit cell. */
  pspatialgeometry sg;
};

/** @brief Representation of parameters. */
typedef struct _parameters parameters;

/** @brief Pointer to a @ref parameters object. */
typedef parameters * pparameters;

/** @brief Range splitting function for the Shifted Potential (SP). 

  @param alpha Damping parameter (void).
  @param real R Cutoff distance.
  @param r Distance between point pair.
  
  @returns Weight of the given radius.
*/
HEADER_PREFIX real 
split_SP (real alpha, real R, real r);

HEADER_PREFIX real
derivative_SP (real alpha, real R, real r);

/** @brief Range splitting function for the Shifted Force potential (SF). 

  @param alpha Damping parameter (void).
  @param real R Cutoff distance.
  @param r Distance between point pair.
  
  @returns Weight of the given radius.*/
HEADER_PREFIX real 
split_SF (real alpha, real R, real r);

HEADER_PREFIX real
derivative_SF (real alpha, real R, real r);

/** @brief Range splitting function for the Damped Shifted Potential (DSP). 

  @param alpha Damping parameter.
  @param real R Cutoff distance.
  @param r Distance between point pair.
  
  @returns Weight of the given radius.*/
HEADER_PREFIX real 
split_DSP (real alpha, real R, real r);

HEADER_PREFIX real
derivative_DSP (real alpha, real R, real r);

/** @brief Range splitting function for the Damped Shifted Force potential (DSF). 

  @param alpha Damping parameter.
  @param real R Cutoff distance.
  @param r Distance between point pair.
  
  @returns Weight of the given radius.*/
HEADER_PREFIX real 
split_DSF (real alpha, real R, real r);

HEADER_PREFIX real
derivative_DSF (real alpha, real R, real r);

/** @brief Range splitting function for the Shifted Potential of second order (SP2). 

  @param alpha Damping parameter (void).
  @param real R Cutoff distance.
  @param r Distance between point pair.
  
  @returns Weight of the given radius.*/
HEADER_PREFIX real 
split_SP2 (real alpha, real R, real r);

HEADER_PREFIX real
derivative_SP2 (real alpha, real R, real r);

/** @brief Range splitting function for the Shifted Potential of third order (SP3). 

  @param alpha Damping parameter (void).
  @param real R Cutoff distance.
  @param r Distance between point pair.
  
  @returns Weight of the given radius.*/
HEADER_PREFIX real 
split_SP3 (real alpha, real R, real r);

HEADER_PREFIX real
derivative_SP3 (real alpha, real R, real r);

/** @brief Coulombic point potential for short range interactions.

@param s Range splitting function.
@param alpha Damping parameter for splitting function.
@param R Cutoff radius.
@param sg Information about the unit cell.
@param x Associated point charge location.
@param y Point charge location where the potential shall be evaluated. 

@returns Value of the short range part of the Coulombic point potential. */
HEADER_PREFIX real 
shortrangePotential_coulomb (split s, real alpha, real R, pspatialgeometry sg, pcreal x, pcreal y);

/** @brief Collection of necessary parameters for coulomb kernel function.

@param s Range splitting function.
@param alpha Damping parameter for splitting function.
@param R Cutoff radius.
@param sg Information about the unit cell.

@returns Parameters for coulomb kernel function. */
HEADER_PREFIX pparameters
setparametersPotential_coulomb (split s, real alpha, real R, pspatialgeometry sg);

/** @brief Kernel function for the short range part of the Coulombic pair potential.

@param x Point charge location.
@param y Point charge location.
@param data Parameters of the Coulomb potential.

@returns Kernel function for Coulombic pair potential, evaluated in x and y. */
HEADER_PREFIX real 
kernel_coulomb (pcreal x, pcreal y, void *data);

/** @brief Gradient of the Coulombic point potential for short range interactions.
 *
 * @param s Range splitting function.
 * @param sd Derivative of the splitting function.
 * @param alpha Damping parameter for splitting function.
 * @param R Cutoff radius.
 * @param sg Information about the unit cell.
 * @param x Associated point charge location.
 * @param y Point charge location where the potential shall be evaluated. 
 * @param f Gradient (with respect to the first variable) of the Coulombic pair 
 *          potential, to be set by this function. */
HEADER_PREFIX void
shortrangeGradient_coulomb (split s, split sd, real alpha, real R, pspatialgeometry sg, 
                            pcreal x, pcreal y, pfield f);

/** @brief Collection of necessary parameters for coulomb gradient function.
 * 
 * @param s Range splitting function.
 * @param sd Derivative of the splitting function.
 * @param alpha Damping parameter for splitting function.
 * @param R Cutoff radius.
 * @param sg Information about the unit cell.
 * 
 * @returns Parameters for coulomb gradient function. */
HEADER_PREFIX pparameters
setparametersGradient_coulomb (split s, split sd, real alpha, real R, pspatialgeometry sg);

/** @brief Gradient function for the short range part of the Coulombic pair potential.
 *
 * @param x Point charge location.
 * @param y Point charge location.
 * @param data Parameters of the Coulomb potential.
 * @param f Gradient (with respect to the first variable) of the Coulombic pair 
 *          potential, to be set by this function. */
HEADER_PREFIX void 
gradient_coulomb (pcreal x, pcreal y, void *data, pfield f);

/** @brief Self-interaction of the point charges.

@param sp Range splitting function.
@param alpha Damping parameter for splitting function.
@param q Vector of point charges in the unit cell.
@param R Cutoff distance.

@returns Value of the self-interaction term of the Coulomb potential. */
HEADER_PREFIX real
selfPotential_coulomb (split sp, real alpha, pavector q, real R);

/** @} */

#endif
