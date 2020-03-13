/** ------------------------------------------------------------
 *    This is the file "tip3p.h" of the KIPS package.
 *    All rights reserved, Jonas Lorenzen 2019
 *  ------------------------------------------------------------ */

/** @file tip3p.h
 *  @author Jonas Lorenzen
 */

#ifndef TIP3P_H
#define TIP3P_H

/** @defgroup tip3p tip3p
 *  @brief Energy and force calculation for the TIP3P model of water.
 *
 *  @attention (Almost) All functions in this module enforce the use of SI units! 
 *            The only exception is the input file for the @ref inputfile_tip3p 
 *            function.
 *  @{ */

#include "basic.h"
#include "spatialgeometry.h"
#include "kernelmatrix.h"
#include "coulomb.h"
#include "lj.h"
#include "rigid.h"

/** @brief Representation of a system of water molecules in the TIP3P model. */
struct _tip3p {
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

/** @brief Type definition of a system of water molecules. */
typedef struct _tip3p tip3p;

/** @brief Pointer to a @ref tip3p object. */
typedef tip3p *ptip3p;


/** ---------------------------------------------------------------------
 *    Global parameters of the TIP3P model
 *  --------------------------------------------------------------------- */

// Zero of the Lennard-Jones potential.
#define SIG_TIP3P 3.15066e-10;

// Depth of the potential well of the Lennard-Jones potential.
#define EPS_TIP3P 1.0566e-21;

// Charge of hydrogen atoms (dimensionless).
#define QH_TIP3P 0.417;

// Oxygen charge (dimensionless).
#define QM_TIP3P -0.834;

// Distance between hydrogen and oxygen atoms.
#define ROH_TIP3P 9.572e-11;

/** Distance between hydrogen atoms, corresponding to an H-O-H angle 
 *  of 104,52°. */
#define RHH_TIP3P 1.5e-10;

/** ---------------------------------------------------------------------
 *    Constructors and destructors
 *  --------------------------------------------------------------------- */
 
/** @brief Create a new empty @ref tip3p object.
 *
 *  Only allocates storage for the components.
 *  Also sets up the vector of charges such that, for each molecule, the 
 *  virtual negative charge is followed by the two positive charges of the 
 *  hydrogen atoms. The site vectors should be ordered accordingly.
 *
 *  @remark Should always be matched by a call to @ref del_tip3p.
 *
 *  @param n Number of water molecules.
 *
 *  @returns Pointer to new tip3p object. */
HEADER_PREFIX ptip3p 
new_tip3p (uint n);

/** @brief Delete a @ref tip3p object.
 *
 *  @param t System of water molecules to be deleted. */
HEADER_PREFIX void
del_tip3p (ptip3p t);
 

/** ---------------------------------------------------------------------
 *    In-/Output
 *  --------------------------------------------------------------------- */

/** @brief Load a system of water molecules from a file. 
 *
 *  Reads the coordinates and types of atoms from a .xyz-type file, 
 *  adds virtual charges corresponding to the TIP3P model and stores 
 *  the system in a @ref tip3p object.
 *
 *  @remark Assumes that the atoms are listed in a way such that the oxygen 
 *  atom of a water molecule is followed by the two corresponding hydrogen 
 *  atoms.
 *  
 *  @attention Assumes the coordinates in the input file to be given in 
 *            units of Angström. Those will then be scaled appropriately.
 *
 *  @param file Name of the input file.
 *
 *  @returns Pointer to a new @ref tip3p object. */
HEADER_PREFIX ptip3p 
inputfile_tip3p (const char *file);

/** ---------------------------------------------------------------------
 *    Energy and force calculation
 *  --------------------------------------------------------------------- */

/** @brief Calculate the potential energy of the TIP3P system.
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
 *  @param sg Geometric information.
 *  @param kc Coulomb potential.
 *  @param klj Lennard-Jones potential.
 *  @param Vc Coulomb potential matrix.
 *  @param Vlj Lennard-Jones potential matrix.
 *  @param t System of water molecules.
 *
 *  @returns Value of the potential energy of the system. */
HEADER_PREFIX real
energy_tip3p (pspatialgeometry sg, pkernelmatrix kc, pkernelmatrix klj, 
                 ph2matrix Vc, ph2matrix Vlj, ptip3p t);

/** @brief Calculate the force vector acting on the center of mass for each 
 *        molecule.
 *
 *  @param fc Concatenated vector of all Coulomb forces.
 *  @param flj Concatenated vector of all Lennard-Jones forces.
 *  @param t System of water molecules. */
HEADER_PREFIX real
calcForce_tip3p (pavector fc, pavector flj, ptip3p t);

/** @brief Calculate the torque vector acting on the center of mass for each 
 *        molecule.
 *
 *  @param fc Concatenated vector of all Coulomb forces.
 *  @param flj Concatenated vector of all Lennard-Jones forces.
 *  @param t System of water molecules. */
HEADER_PREFIX real
calcTorque_tip3p (pavector fc, pavector flj, ptip3p t);

/** @brief Recalculation of positions, forces and torques.
 *
 *  Computes the new absolute atomic positions, when the center of mass positions 
 *  and relative atomic positions have been altered. Then updates forces and 
 *  torques for each molecule. 
 *
 *  @param sg Geometric information.
 *  @param kc Coulomb potential.
 *  @param klj Lennard-Jones potential.
 *  @param Fc Coulomb gradient matrix.
 *  @param Flj Lennard-Jones gradient matrix.
 *  @param fc Concatenated vector of all Coulomb forces.
 *  @param flj Concatenated vector of all Lennard-Jones forces.
 *  @param t System of water molecules. */
HEADER_PREFIX void
update_tip3p (pspatialgeometry sg, pkernelmatrix kc, pkernelmatrix klj, 
              ph2matrix Fc, ph2matrix Flj, pavector fc, pavector flj, ptip3p t);

/** ---------------------------------------------------------------------
 *    Optimization
 *  --------------------------------------------------------------------- */

/** @brief Gradient descent method for minimization of the potential energy.
 *
 *  Does the same calculation as @ref energy_tip3p, then performs one step of 
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
 *  @param Vc Coulomb potential matrix.
 *  @param Vlj Lennard-Jones potential matrix.
 *  @param Fc Coulomb gradient matrix.
 *  @param Flj Lennard-Jones gradient matrix.
 *  @param lambda Initial step size for the line search.
 *  @param nu Factor by which the step size is decreased during the line search.
 *  @param m Maximal number of line search steps.
 *  @param t System of water molecules.
 *
 *  @returns New value of the potential energy of the system. */
HEADER_PREFIX real
gradientDescent_tip3p (pspatialgeometry sg, pkernelmatrix kc, pkernelmatrix klj, 
                       ph2matrix Vc, ph2matrix Fc, ph2matrix Vlj, ph2matrix Flj, 
                       real lambda, real nu, uint m, ptip3p t);

//HEADER_PREFIX real
//gradientDescentRotation_tip3p (pspatialgeometry sg, pkernelmatrix kc, pkernelmatrix klj, 
//                              ph2matrix Vc, ph2matrix Fc, ph2matrix Vlj, ph2matrix Flj, 
//                              real lambda, real nu, uint m, ptip3p t);

/** ----------------------------------------------------------------------
 *    Time-stepping methods
 *  ---------------------------------------------------------------------- */

/** @brief Construct the inertia matrix and calculate reference orientation for each 
 *        molecule.
 *
 *  @param t System of water molecules. */
HEADER_PREFIX void
initMolecules_tip3p (ptip3p t);

/** @brief One step of the Velocity Verlet scheme.
 *
 *  Performs translation and reorientation of the rigid molecules. Then recalculates 
 *  forces and torques and finally updates velocities of the translational and 
 *  rotational motion, respectively.
 *
 *  @param sg Geometric information.
 *  @param kc Coulomb potential.
 *  @param klj Lennard-Jones potential.
 *  @param Fc Coulomb gradient matrix.
 *  @param Flj Lennard-Jones gradient matrix.
 *  @param fc Concatenated vector of all Coulomb forces.
 *  @param flj Concatenated vector of all Lennard-Jones forces.
 *  @param delta Timestep length.
 *  @param t System of water molecules. 
 *
 *  @returns Current kinetic energy of the system. */
HEADER_PREFIX real
velocityVerlet_tip3p (pspatialgeometry sg, pkernelmatrix kc, pkernelmatrix klj, 
                      ph2matrix Fc, ph2matrix Flj, pavector fc, pavector flj, 
                      real delta, ptip3p t);

/** @} */

#endif
