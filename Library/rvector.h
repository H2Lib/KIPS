
/* ------------------------------------------------------------
 * This is the file "rvector.h" of the KIPS package.
 * All rights reserved, Steffen Boerm 2015
 * ------------------------------------------------------------ */

/** @file rvector.h
 *  @author Steffen B&ouml;rm
 */

#ifndef RVECTOR_H
#define RVECTOR_H

/** @defgroup rvector rvector
 *  @brief Representation of a real-valued vector as an array.
 *
 *  The @ref rvector class is used to store singular values and
 *  eigenvalues of self-adjoint matrices.
 *  @{ */

/** Representation of a vector as an array. */
typedef struct _rvector rvector;

/** Pointer to a @ref rvector object. */
typedef rvector *prvector;

/** Pointer to a constant @ref rvector object. */
typedef const rvector *pcrvector;

#include <assert.h>
#include <math.h>
#include <stdlib.h>

#include "settings.h"

/** Representation of a real-valued vector as an array. */
struct _rvector {
  /** @brief Vector coefficients. */
  real *v;

  /** @brief Vector dimension. */
  uint size;
};

/* ------------------------------------------------------------
   Constructors and destructors
   ------------------------------------------------------------ */

/** @brief Initialize an @ref rvector object.
 *
 *  Sets up the components of the object and allocates storage for
 *  the coefficient array.
 *
 *  @remark Should always be matched by a call to @ref uninit_rvector.
 *
 *  @param v Object to be initialized.
 *  @param size Dimension of the new vector.
 *  @returns Initialized @ref avector object. */
HEADER_PREFIX prvector
init_rvector(prvector v, uint size);

/** @brief Uninitialize an @ref rvector object.
 *
 *  Invalidates pointers, freeing corresponding storage,
 *  and prepares the object for deletion.
 *
 *  @param v Object to be uninitialized. */
HEADER_PREFIX void
uninit_rvector(prvector v);

/** @brief Create a new @ref rvector object.
 *
 *  Allocates storage for the object an sets up its components.
 *
 *  @remark Should always be matched by a call to @ref del_rvector.
 *
 *  @param size Dimension of the new vector.
 *  @returns New @ref avector object. */
HEADER_PREFIX prvector
new_rvector(uint size);

/** @brief Delete an @ref rvector object.
 *
 *  Releases the storage corresponding to the object.
 *
 *  @param v Object to be deleted. */
HEADER_PREFIX void
del_rvector(prvector v);

/** @brief Change the dimension of an @ref rvector object without
 *  preserving its coefficients.
 *
 *  Allocates new storage for the coefficients and releases the
 *  old storage.
 *
 *  @param v Vector to be resized.
 *  @param size New vector dimension. */
HEADER_PREFIX void
resize_rvector(prvector v, uint size);

/* ------------------------------------------------------------
   Access methods
   ------------------------------------------------------------ */

#ifdef __GNUC__
INLINE_PREFIX real
getentry_rvector(pcrvector, uint) __attribute__((const,unused));
INLINE_PREFIX void
setentry_rvector(prvector, uint, real) __attribute__((unused));
#endif

/** @brief Read a vector entry @f$v_i@f$.
 *
 *  @param v Vector @f$v@f$.
 *  @param i Index @f$i@f$.
 *  @returns Vector entry @f$v_i@f$. */
INLINE_PREFIX real
getentry_rvector(pcrvector v, uint i)
{
#ifdef FULL_DEBUG
  assert(i < v->dim);
#endif
  return v->v[i];
}

/** @brief Set a vector entry, @f$v_i \gets x@f$.
 *
 *  @param v Vector @f$v@f$.
 *  @param i Index @f$i@f$.
 *  @param x New value of @f$v_i@f$. */
INLINE_PREFIX void
setentry_rvector(prvector v, uint i, real x)
{
#ifdef FULL_DEBUG
  assert(i < v->dim);
#endif

  v->v[i] = x;
}

/* ------------------------------------------------------------
   Statistics
   ------------------------------------------------------------ */

/** @brief Get number of currently initialized @ref rvector objects.
 *
 *  Calls to initialization functions like @ref init_rvector and
 *  constructors like @ref new_rvector increase an internal counter,
 *  while @ref uninit_rvector and @ref del_rvector decrease it.
 *
 *  @remark Use this function to check whether a program correctly cleans
 *  up temporary variables.
 *
 *  @returns Number of currently initialized @ref rvector objects. */
HEADER_PREFIX uint
getactives_rvector();

/** @brief Get size of a given @ref rvector object.
 *
 *  Computes the size of the @ref rvector object and the storage
 *  allocated for the coefficients.
 *
 *  @param v Vector object.
 *  @returns Size of allocated storage in bytes. */
HEADER_PREFIX size_t
getsize_rvector(pcrvector v);

/* ------------------------------------------------------------
   Simple utility functions
   ------------------------------------------------------------ */

/** @brief Copy a vector into another vector, @f$w \gets v@f$.
 *
 *  The sizes of both vectors have to be identical.
 *
 *  @param v Source vector.
 *  @param w Target vector. */
HEADER_PREFIX void
copy_rvector(pcrvector v, prvector w);

/** @brief Print a vector.
 *
 *  @param v Vector @f$v@f$. */
HEADER_PREFIX void
print_rvector(pcrvector v);

/** @} */

#endif
