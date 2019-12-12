/* ------------------------------------------------------------
 * This is the file "quaternion.h" of the KIPS package.
 * All rights reserved, Jonas Lorenzen 2019
 * ------------------------------------------------------------ */

/** @file quaternion.h
 *  @author Jonas Lorenzen
 */

#ifndef QUATERNION_H
#define QUATERNION_H

/** @defgroup quaternion quaternion
 *  @brief 3D-rotation and matrix diagonalization using quaternions.
 *  @{ */

#include "math.h"
#include "amatrix.h"

/** @brief Pointer to a quaternion. (Identical to a real pointer that 
 *        points to exactly four real variables.) */
typedef preal pquaternion;

/** @brief Pointer to a constant quaternion. */
typedef pcreal pcquaternion;


/* ---------------------------------------------------------------------
 *    Constructors and destructors
 * --------------------------------------------------------------------- */

/** @brief Generate a new @ref quaternion object.
 *
 *  @returns Pointer to new quaternion object. */
HEADER_PREFIX pquaternion 
new_quaternion ();

/** @brief Delete a @ref quaternion object. 
 *
 *  @param q Quaternion to be deleted. */
HEADER_PREFIX void
del_quaternion (pquaternion q);

/* ---------------------------------------------------------------------
 *    Basic operations
 * --------------------------------------------------------------------- */

/** @brief Normalize a quaternion.
 *
 *  Divide the quaternion components by the euclidean norm of the 
 *  quaternion in order to retain a unit quaternion.
 *
 *  @param q Quaternion that has to be normalized. */
HEADER_PREFIX void
normalize_quaternion (pquaternion q);

/** @brief Quaternion right multiplication.
 *
 *  Performs the operation @f$q \gets qp@f$.
 *
 *  @param p Multiplier.
 *  @param q In- and output quaternion. */
HEADER_PREFIX void 
update_quaternion (pcquaternion p, pquaternion q);

/** @brief Build a 3D-rotation matrix based on a quaternion.
 *
 *  Fills a 3x3-matrix with the correct entries. 
 *  Does not allocate (or free) memory.
 *
 *  @param q Rotation quaternion.
 *  @param R Rotation matrix. */
HEADER_PREFIX void 
buildRotation_quaternion (pcquaternion q, pamatrix R);


/* ---------------------------------------------------------------------
 *    Diagonalization
 * --------------------------------------------------------------------- */

/** @brief Jacobi eigenvalue iteration with quaternions.
 *
 *  Performs the Jacobi diagonalization of a self-adjoint real 
 *  or complex 3x3-matrix. The sequence of Givens rotations used 
 *  are represented by a unit quaternion.
 *
 *  @param q Rotation quaternion.
 *  @param R Rotation matrix. */
HEADER_PREFIX uint 
Jacobi_quaternion (pamatrix A, pquaternion q);

#endif
