/* ------------------------------------------------------------
 * This is the file "tip4p.h" of the KIPS package.
 * All rights reserved, Jonas Lorenzen 2019
 * ------------------------------------------------------------ */

/** @file tip4p.h
 *  @author Jonas Lorenzen
 */

#ifndef TIP4P_H
#define TIP4P_H

/** @defgroup tip4p tip4p
 *  @brief Energy and force calculation for the TIP4P model of water.
 *  @{ */

#include "basic.h"
#include "spatialgeometry.h"
#include "kernelmatrix.h"
#include "coulomb.h"
#include "lj.h"
#include "rigid.h"

/** @brief Representation of a system of water molecules in the TIP4P model. */
struct _tip4p {
  /** @brief Number of water molecules. */
  uint n;
  
  /** @brief Array of molecules. */
  pmolecule *mol;
  
  /** @brief Array of charged sites. */
  preal *xc;
  
  /** @brief Molecule numbers for charged sites from 1 to n. */
  uint *xcmol;
  
  /** @brief Vector of charges. */
  pavector q;
  
  /** @brief Array of uncharged sites. */
  preal *xn;
  
  /** @brief Molecule numbers for uncharged sites from 1 to n. */
  uint *xnmol;
};

/** @brief Representation of a system of water molecules. */
typedef struct _tip4p tip4p;

/** @brief Pointer to a @ref tip4p object. */
typedef tip4p *ptip4p;


/* ---------------------------------------------------------------------
 *    Global parameters of the TIP4P model
 * --------------------------------------------------------------------- */

/** Lengths are measured in Angström. 
 *  Charges are in terms of the elementary charge.*/

// Zero of the Lennard-Jones potential.
#define SIG_TIP4P 3.154;

// Depth of the potential well in Joule.
#define EPS_TIP4P 1.0769e-21;

// Charge of hydrogen atoms.
#define QH_TIP4P 0.52;

// Virtual charges.
#define QM_TIP4P -1.04;

// Distance between hydrogen and oxygen atoms.
#define ROH_TIP4P 0.9572;

// Distance between oxygen atoms and virtual charges.
#define ROM_TIP4P 0.15;

/** Distance between hydrogen atoms, corresponding to an H-O-H angle 
 *  of 104,52°. */
#define RHH_TIP4P 1.5;

/* ---------------------------------------------------------------------
 *    Constructors and destructors
 * --------------------------------------------------------------------- */
 
/** @brief Create a new empty @ref tip4p object.
 *
 *  Only allocates storage for the components.
 *  Also sets up the vector of charges such that, for each molecule, the 
 *  virtual negative charge is followed by the two positive charges of the 
 *  hydrogen atoms. The site vectors should be ordered accordingly.
 *
 *  @remark Should always be matched by a call to @ref del_tip4p.
 *
 *  @param n Number of water molecules.
 *
 *  @returns Pointer to new tip4p object. */
HEADER_PREFIX ptip4p 
new_tip4p (uint n);

/** @brief Delete a @ref tip4p object.
 *
 *  Releases the storage corresponding to the given tip4p object.
 *
 *  @param t Object to be deleted. */
HEADER_PREFIX void
del_tip4p (ptip4p t);
 

/* ---------------------------------------------------------------------
 *    In-/Output
 * --------------------------------------------------------------------- */

/** @brief Load a system of water molecules from a file. 
 *
 *  Reads the coordinates and types of atoms from a .xyz-type file, 
 *  adds virtual charges corresponding to the TIP4P model and stores 
 *  the system in a @ref tip4p object.
 *
 *  @remark: Assumes that the atoms are listed in a way such that the oxygen 
 *  atom of a water molecule is followed by the two corresponding hydrogen 
 *  atoms.
 *  Furthermore assumes all distances to be measured in Angström.
 *
 *  @param file Name of the input file.
 *
 *  @returns Pointer to a new tip4p object corresponding to the input file. */
HEADER_PREFIX ptip4p 
inputfile_tip4p (const char *file);

HEADER_PREFIX void 
calc_virtual_charge (pcreal h1, pcreal h2, pcreal o, preal m);

/* ---------------------------------------------------------------------
 *    Energy and force calculation
 * --------------------------------------------------------------------- */

/** @brief Calculate the potential energy of the TIP4P system.
 *
 *  Calculate the electrostatic energy via smooth cutoff schemes for the 
 *  Coulomb and Lennard-Jones potentials, respectively. 
 *  An @f$\mathcal{H}^2@f$-matrix structure with periodic boundary conditions 
 *  is used for each of the two potentials. The raw geometric and algebraic 
 *  structure should already be given by the corresponding @ref spatialgeometry, 
 *  @ref kernelmatrix and @ref h2matrix objects, respectively. 
 *  This function will then fill the matrices with the suitable values of the 
 *  potentials and calculate the potential energy as a sum over all entries.
 *
 *  @remark All distances are supposed to be measured in Angström.
 *
 *  @param sg Geometric information.
 *  @param kc Coulomb potential.
 *  @param klj Lennard-Jones potential.
 *  @param Vc @f$\mathcal{H}^2@f$-matrix for the Coulomb potential.
 *  @param Vlj @f$\mathcal{H}^2@f$-matrix for the Lennard-Jones potential.
 *  @param t System of water molecules.
 *
 *  @returns Value of the potential energy of the system in J/mol. */
HEADER_PREFIX real
energy_tip4p (pspatialgeometry sg, pkernelmatrix kc, pkernelmatrix klj, 
                 ph2matrix Vc, ph2matrix Vlj, ptip4p t);

/** @brief Calculate the force vector acting on the center of mass for each 
 *        molecule.
 *
 *  @remark All distances are supposed to be measured in Angström.
 *
 *  @param fc Concatenated vector of all Coulomb forces.
 *  @param flj Concatenated vector of all Lennard-Jones forces.
 *  @param t System of water molecules. */
HEADER_PREFIX real
calcForce_tip4p (pavector fc, pavector flj, ptip4p t);

/** @brief Calculate the torque vector acting on the center of mass for each 
 *        molecule.
 *
 *  @remark All distances are supposed to be measured in Angström.
 *
 *  @param fc Concatenated vector of all Coulomb forces.
 *  @param flj Concatenated vector of all Lennard-Jones forces.
 *  @param t System of water molecules. */
HEADER_PREFIX real
calcTorque_tip4p (pavector fc, pavector flj, ptip4p t);


/* ---------------------------------------------------------------------
 *    Optimization
 * --------------------------------------------------------------------- */

/** @brief Gradient descent method for minimization of the potential energy.
 *
 *  Does the same calculation as @ref energy_tip4p, then performs one step of 
 *  the gradient descent method and finally recalculates the potential energy.
 *  The step size is determined by a backtracking line search. A control 
 *  parameter of one half for the decrease of energy is used.
 *
 *  @remark Currently only optimizes the translational component. 
 *          Rotation is yet to be implemented.
 *
 *  @param sg Geometric information.
 *  @param kc Coulomb potential.
 *  @param klj Lennard-Jones potential.
 *  @param Vc @f$\mathcal{H}^2@f$-matrix for the Coulomb potential.
 *  @param Vlj @f$\mathcal{H}^2@f$-matrix for the Lennard-Jones potential.
 *  @param Fc Gradient matrix of the Coulomb potential.
 *  @param Flj Gradient matrix of the Lennard-Jones potential.
 *  @param lambda Initial step size for the line search.
 *  @param nu Factor by which the step size is decreased during the line search.
 *  @param m Maximal number of line search steps.
 *  @param t System of water molecules.
 *
 *  @returns New value of the potential energy of the system in J/mol. */
HEADER_PREFIX real
gradientDescent_tip4p (pspatialgeometry sg, pkernelmatrix kc, pkernelmatrix klj, 
                       ph2matrix Vc, ph2matrix Fc, ph2matrix Vlj, ph2matrix Flj, 
                       real lambda, real nu, uint m, ptip4p t);

HEADER_PREFIX real
gradientDescentRotation_tip4p (pspatialgeometry sg, pkernelmatrix kc, pkernelmatrix klj, 
                              ph2matrix Vc, ph2matrix Fc, ph2matrix Vlj, ph2matrix Flj, 
                              real lambda, real nu, uint m, ptip4p t);

/* ----------------------------------------------------------------------
 *    Rotational motion
 * ---------------------------------------------------------------------- */

/** @brief Construct the inertia matrix for each molecule.
 *
 *  @remark Assumes lengths to be measured in Angström. 
 *          The diagonalization, e. g. via @ref Jacobi_quaternion, is invariant 
 *          to scaling anyway, but the principal moments of inertia might be scaled 
 *          differently in other units.
 *
 *  @param t System of water molecules. */
HEADER_PREFIX void
inertia_tip4p (ptip4p t);

/* ----------------------------------------------------------------------
 *    Translational motion
 * ---------------------------------------------------------------------- */

#endif
