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
 *  This is a weight function for point distances up to a certain Cutoff radius.*/
typedef real (*split) (real alpha, real R, real r);

/** @brief Representation of the necessary parameters for a kernel function 
 *  that evaluates a coulombic potential or forces. */
struct _coulomb {
  /** @brief Range splitting function. */
  split s;
  
  /** @brief Derivative of the splitting function. */
  split sd;
  
  /** @brief Damping parameter for certain splitting functions. */
  real alpha;
  
  /** @brief Cutoff distance. */
  real R;
  
  /** @brief Information about the unit cell. */
  pcspatialgeometry sg;
};

/** @brief Representation of parameters. */
typedef struct _coulomb coulomb;

/** @brief Pointer to a parameter object. */
typedef coulomb * pcoulomb;


/*  ----------------------------------------------------------------
 *    Range splitting functions for cutoff schemes
 *  ---------------------------------------------------------------- */

/** @brief Range splitting function for the simple cutoff potential. 
 *
 *  @param alpha Damping parameter (void).
 *  @param real R Cutoff distance.
 *  @param r Distance between point pair.
 * 
 *  @returns Weight of the given radius. */
HEADER_PREFIX real 
split_cutoff (real alpha, real R, real r);

/** @brief Derivative of the range splitting function for the simple 
 *        cutoff potential. 
 *
 *  @param alpha Damping parameter (void).
 *  @param real R Cutoff distance.
 *  @param r Distance between point pair.
 * 
 *  @returns Weight of the given radius. */
HEADER_PREFIX real
derivative_cutoff (real alpha, real R, real r);

/** @brief Range splitting function for the Shifted Potential (SP). 
 *
 *  @param alpha Damping parameter (void).
 *  @param real R Cutoff distance.
 *  @param r Distance between point pair.
 * 
 *  @returns Weight of the given radius. */
HEADER_PREFIX real 
split_SP (real alpha, real R, real r);

/** @brief Derivative of the range splitting function for the Shifted 
 *        Potential (SP). 
 *
 *  @param alpha Damping parameter (void).
 *  @param real R Cutoff distance.
 *  @param r Distance between point pair.
 * 
 *  @returns Weight of the given radius. */
HEADER_PREFIX real
derivative_SP (real alpha, real R, real r);

/** @brief Range splitting function for the Shifted Force potential (SF). 
 *  
 *  @param alpha Damping parameter (void).
 *  @param real R Cutoff distance.
 *  @param r Distance between point pair.
 *
 *  @returns Weight of the given radius. */
HEADER_PREFIX real 
split_SF (real alpha, real R, real r);

/** @brief Derivative of the range splitting function for the Shifted 
 *        Force potential (SF). 
 *  
 *  @param alpha Damping parameter (void).
 *  @param real R Cutoff distance.
 *  @param r Distance between point pair.
 *
 *  @returns Weight of the given radius. */
HEADER_PREFIX real
derivative_SF (real alpha, real R, real r);

/** @brief Range splitting function for the Damped Shifted Potential (DSP). 
 *
 *  @param alpha Damping parameter.
 *  @param real R Cutoff distance.
 *  @param r Distance between point pair.
 * 
 *  @returns Weight of the given radius. */
HEADER_PREFIX real 
split_DSP (real alpha, real R, real r);

/** @brief Derivative of the range splitting function for the Damped 
 *        Shifted Potential (DSP). 
 *
 *  @param alpha Damping parameter.
 *  @param real R Cutoff distance.
 *  @param r Distance between point pair.
 * 
 *  @returns Weight of the given radius. */
HEADER_PREFIX real
derivative_DSP (real alpha, real R, real r);

/** @brief Range splitting function for the Damped Shifted Force 
 *        potential (DSF). 
 *
 *  @param alpha Damping parameter.
 *  @param real R Cutoff distance.
 *  @param r Distance between point pair.
 * 
 *  @returns Weight of the given radius. */
HEADER_PREFIX real 
split_DSF (real alpha, real R, real r);

/** @brief Derivative of the range splitting function for the Damped 
 *        Shifted Force potential (DSF). 
 *
 *  @param alpha Damping parameter.
 *  @param real R Cutoff distance.
 *  @param r Distance between point pair.
 * 
 *  @returns Weight of the given radius. */
HEADER_PREFIX real
derivative_DSF (real alpha, real R, real r);

/** @brief Range splitting function for the Shifted Potential of 
 *        second order (SP2). 
 *
 *  @param alpha Damping parameter (void).
 *  @param real R Cutoff distance.
 *  @param r Distance between point pair.
 * 
 *  @returns Weight of the given radius. */
HEADER_PREFIX real 
split_SP2 (real alpha, real R, real r);

/** @brief Derivative of the range splitting function for the Shifted 
 *        Potential of second order (SP2). 
 *
 *  @param alpha Damping parameter (void).
 *  @param real R Cutoff distance.
 *  @param r Distance between point pair.
 * 
 *  @returns Weight of the given radius. */
HEADER_PREFIX real
derivative_SP2 (real alpha, real R, real r);

/** @brief Range splitting function for the Shifted Potential of 
 *        third order (SP3). 
 *
 *  @param alpha Damping parameter (void).
 *  @param real R Cutoff distance.
 *  @param r Distance between point pair.
 * 
 *  @returns Weight of the given radius.*/
HEADER_PREFIX real 
split_SP3 (real alpha, real R, real r);

/** @brief Derivative of the range splitting function for the 
 *        Shifted Potential of third order (SP3). 
 *
 *  @param alpha Damping parameter (void).
 *  @param real R Cutoff distance.
 *  @param r Distance between point pair.
 * 
 *  @returns Weight of the given radius.*/
HEADER_PREFIX real
derivative_SP3 (real alpha, real R, real r);


/*  ----------------------------------------------------------------
 *    Geometric data
 *  ---------------------------------------------------------------- */

/** @brief Minimal diameter for @ref spatialcluster objects in order to 
 *  ensure that intramolecular interactions are ignored in the calculation 
 *  of Coulomb forces.
 *
 *  @param eta Admissibility parameter for the corresponding @ref h2matrix 
 *        object.
 *  @param rmol Length of the largest molecule considered in the simulation.
 *
 *  @returns Minimal diameter needed for construction of a spatial cluster tree. */
HEADER_PREFIX real
mindiam_coulomb (real eta, real rmol);


/*  ----------------------------------------------------------------
 *    Coulomb Potential
 *  ---------------------------------------------------------------- */

/** @brief Collection of necessary parameters for coulomb kernel function.
 *
 *  @param s Range splitting function.
 *  @param alpha Damping parameter for splitting function.
 *  @param R Cutoff radius.
 *  @param sg Information about the unit cell.
 *
 *  @returns Parameters for coulomb kernel function. */
HEADER_PREFIX pcoulomb
setparametersPotential_coulomb (split s, real alpha, real R, pcspatialgeometry sg);

/** @brief Kernel function for the short-range part of the Coulombic pair potential.
 *
 *  @remark The calculation is dimensionless and free of constants such as 
 *          vacuum permittivity. The considered charge of particles is +1. This can 
 *          be modified by multiplication with a suitable charge vector.
 *
 *  @param x Evaluation point.
 *  @param y Evaluation point.
 *  @param xmol Molecule number of x.
 *  @param ymol Molecule number of y.
 *  @param data Parameters of the Coulomb potential.
 *
 *  @returns Kernel function for Coulombic pair potential, evaluated in x and y. */
HEADER_PREFIX real 
kernel_coulomb (pcreal x, pcreal y, uint xmol, uint ymol, void *data);

/** @brief Self-term of the potential energy.
 *
 *  This term counteracts self-interactions of image charges due to smooth 
 *  cutoff schemes.
 *
 *  @param sp Range splitting function.
 *  @param alpha Damping parameter for splitting function.
 *  @param q Vector of point charges in the unit cell.
 *  @param R Cutoff distance.
 *
 *  @returns Value of the self-term of the Coulomb potential. */
HEADER_PREFIX real
selfterm_coulomb (split sp, real alpha, pcavector q, real R);

/** @brief Potential energy of a given system of charged sites.
 *
 *  @remark The input potential matrix is supposed to have the correct 
 *          structure, but leaf and nearfield matrices will be updated 
 *          within this function.
 *          Assumes all distances to be measured in Angström.
 *
 *  @param q (Dimensionless) Charge vector.
 *  @param kc Kernel matrix of the Coulomb potential.
 *  @param sg Geometric information.
 *  @param Vc Potential matrix of @ref h2matrix type.
 *
 *  @returns Value of the potential energy in units of Joule/mol. */
HEADER_PREFIX real
energy_coulomb (pcavector q, pckernelmatrix kc, pspatialgeometry sg, ph2matrix Vc);


/*  ----------------------------------------------------------------
 *    Coulomb Forces
 *  ---------------------------------------------------------------- */

/** @brief Collection of necessary parameters for coulomb gradient function.
 * 
 *  @param s Range splitting function.
 *  @param sd Derivative of the splitting function.
 *  @param alpha Damping parameter for splitting function.
 *  @param R Cutoff radius.
 *  @param sg Information about the unit cell.
 * 
 *  @returns Parameters for coulomb gradient function. */
HEADER_PREFIX pcoulomb
setparametersGradient_coulomb (split s, split sd, real alpha, real R, 
                               pcspatialgeometry sg);

/** @brief Gradient function for the short-range part of the Coulombic pair potential.
 *
 *  @remark The calculation is dimensionless and free of constants such as 
 *          vacuum permittivity. The considered charge of particles is +1. This can 
 *          be modified by multiplication with a suitable charge vector.
 *
 *  @param x Point charge location.
 *  @param y Point that exerts force.
 *  @param xmol Molecule number of x.
 *  @param ymol Molecule number of y.
 *  @param data Parameters of the Coulomb potential.
 *  @param f Gradient of the Coulombic pair potential, to be set by this function. */
HEADER_PREFIX void 
gradient_coulomb (pcreal x, pcreal y, uint xmol, uint ymol, void *data, pfield f);

/** @brief Compute a concatenated vector of the forces acting on the charged sites 
 *        of the given system.
 *
 *  The i-th to (i+dim-1)-th entries of the vector correspond to the force vector 
 *  acting on the i-th point charge.
 *
 *  @remark The input gradient matrix is supposed to have the correct 
 *          structure, but leaf and nearfield matrices will be updated 
 *          within this function.
 *          Assumes all distances to be measured in Angström.
 *
 *  @param q Vector of charges.
 *  @param kc Kernel matrix for the Coulomb potential.
 *  @param sg Geometric information.
 *  @param Fc Gradient matrix of @ref h2matrix type. 
 *  @param f Pointer to force vector. */
HEADER_PREFIX void
fullForce_coulomb (pcavector q, pckernelmatrix kc, pspatialgeometry sg, 
                   ph2matrix Fc, pavector f);

/** @} */

#endif
