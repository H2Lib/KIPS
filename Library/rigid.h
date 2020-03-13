/** ------------------------------------------------------------
 *    This is the file "rigid.h" of the KIPS package.
 *    All rights reserved, Jonas Lorenzen 2019
 *  ------------------------------------------------------------ */

/** @file rigid.h
 *  @author Jonas Lorenzen
 */

#ifndef RIGID_H
#define RIGID_H

/** @defgroup rigid rigid
 *  @brief Methods for solution of the Newton-Euler equations of motion of 
 *        rigid bodies in three-dimensional space.
 *  @{ */

#include "basic.h"
#include "quaternion.h"
#include "amatrix.h"
#include "avector.h"

/** @brief Representation of a rigid molecule. */
struct _molecule {
  /** @brief Number of atoms. */
  uint n;
  
  /** @brief Mass of the molecule. */
  real m;
  
  /** @brief Coordinates of the center of mass. */
  preal com;
  
  /** @brief Relative coordinates of the atoms. */
  preal *r;
  
  /** @brief Relative coordinates in the reference orientation. */
  preal *rp;
  
  /** @brief Force vector acting on the center of mass. */
  preal f;
  
  /** @brief Translational velocity vector. */
  preal v;
  
  /** @brief Torque vector acting on the center of mass. */
  preal t;
  
  /** @brief Torque vector in the reference orientation. */
  preal tp;
  
  /** @brief Angular velocity vector in the reference frame. */
  preal w;
  
  /** @brief Memory that may be used initially for the inertia matrix 
   *        and later on for a rotation matrix. 
   *        Matrices will be stored compatible with BLAS 
   *        (in column-major order). */
  preal R;
  
  /** @brief Principal moments of inertia. */
  preal I;
  
  /** @brief Quaternion that describes the molecule's orientation in 
   *        three-dimensional space. */
  pquaternion q;
  
  /** @brief Quaternion velocity. */
  pquaternion qv;
  
  /** @brief Quaternion acceleration. */
  pquaternion qa;
  
  /** @brief Intermediary storage for half-step quantities. */
  pquaternion qh;
};

/** @brief Type definition of a rigid molecule. */
typedef struct _molecule molecule;

/** @brief Pointer to a @ref molecule object. */
typedef molecule *pmolecule;


/** ---------------------------------------------------------------------
 *    Constructors and destructors
 *  --------------------------------------------------------------------- */

/** @brief Create a new empty @ref molecule object.
 *
 *  Only allocates storage for the components.
 *  
 *  @remark Should always be matched by a call to @ref del_molecule.
 *
 *  @param n Number of atoms.
 *
 *  @returns Pointer to new @ref molecule object. */
HEADER_PREFIX pmolecule
new_molecule (uint n);

/** @brief Delete a @ref molecule object.
 *
 *  @param mol Molecule to be deleted. */
HEADER_PREFIX void
del_molecule (pmolecule mol);


/** ---------------------------------------------------------------------
 *    Reference orientation
 *  --------------------------------------------------------------------- */

HEADER_PREFIX void
initRefOrientationFromFile_molecule (const char *file, uint n, pmolecule *mol);

/** @brief Compute reference orientation of a rigid molecule.
 *
 *  Diagonalizes the inertia matrix to obtain principal moments of 
 *  inertia, relative coordinates of the reference orientation and a 
 *  quaternion corresponding to the suitable rotation as well as the 
 *  corresponding rotation matrix. Furthermore transforms the relative 
 *  atom positions and the torque towards the reference orientation.
 *  Uses the Jacobi eigenvalue algorithm.
 *
 *  @remark Assumes that the matrix slot of the @ref molecule object 
 *          is initialized with the correct inertia matrix of the 
 *          molecule. This will be overwritten by the rotation 
 *          matrix corresponding to the reference orientation.
 *
 *  @param mol Underlying molecule. */
HEADER_PREFIX void
refOrientation_molecule (pmolecule mol);

/** @brief Compute the new orientation based on a quaternion rotation.
 *
 *  Calculates new rotation matrix based on a given quaternion and 
 *  rotates the molecule accordingly.
 *
 *  @param q Quaternion defining the rotation.
 *  @param mol Underlying molecule. */
HEADER_PREFIX void
adjust_molecule (pcquaternion q, pmolecule mol);


/** ---------------------------------------------------------------------
 *    Time-stepping methods
 *  --------------------------------------------------------------------- */

HEADER_PREFIX real 
initVelocitiesFromFile_molecule (const char *file, uint n, pmolecule *mol);

/** @brief Perform the position integration of the Velocity Verlet method.
 *
 *  Calculates the new position of the center of mass and the new 
 *  orientation of the molecule after a given timestep.
 *
 *  @param delta Length of the timestep.
 *  @param mol Underlying molecule. */
HEADER_PREFIX void
positionVerlet_molecule (real delta, pmolecule mol);

/** @brief Perform the velocity integration of the Velocity Verlet method.
 *
 *  Calculates the new velocities of the translational and rotational motion, 
 *  respectively, after a given timestep.
 *
 *  @param delta Length of the timestep.
 *  @param mol Underlying molecule. 
 *
 *  @returns Current kinetic energy of the molecule. */
HEADER_PREFIX real
velocityVerlet_molecule (real delta, pmolecule mol);

/** @} */

#endif
