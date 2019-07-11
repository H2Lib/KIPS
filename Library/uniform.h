
/* ------------------------------------------------------------
 * This is the file "uniform.h" of the KIPS package.
 * All rights reserved, Steffen Boerm 2009
 * ------------------------------------------------------------ */

/** @file uniform.h
 *  @author Steffen B&ouml;rm */

#ifndef UNIFORM_H
#define UNIFORM_H

/** \defgroup uniform uniform
 *  @brief Representation of an admissible block for
 *  @f$ \mathcal H^2 @f$-matrices.
 *
 *  In case of @f$ \mathcal H^2@f$-matrices an admissible subblock is
 *  represented by
 *  @f[
 *  G|_{t \times s} = V_t S_b W_s^*,
 *  @f]
 *  where @f$ V_t, W_s @f$ are the corresponding clusterbasis.
 *
 *  A @ref _uniform "uniform" object stores the coupling-matrix @f$ S_b @f$
 *  as well as references to the row- and column-clusterbasis.
 *  @{ */

/** @brief Representation of an admissible block for
 *  @f$ \mathcal H^2 @f$-matrices. */
typedef struct _uniform uniform;

/** @brief Pointer to @ref _uniform "uniform" object. */
typedef uniform *puniform;

/** @brief Pointer to a constant @ref _uniform "uniform" object. */
typedef const uniform *pcuniform;

#include "settings.h"
#include "amatrix.h"
#include "clusterbasis.h"

/** @brief Representation of an admissible block for
 *  @f$ \mathcal H^2 @f$-matrices.
 *
 *  <tt>rb</tt> and <tt>cb</tt> are references to the row- and column-
 *  @ref _clusterbasis "clusterbasis" respectively.
 *  The coupling-matrix is stored inside an @ref _amatrix "amatrix" <tt>S</tt>.
 *
 *  For convenience reasons there a two lists storing all @ref _uniform "uniform"
 *  objects belonging to the same block row defined by the row @ref _clusterbasis
 *  "clusterbasis" <tt>rb</tt> via <tt>rnext, rprev</tt> or to the same block
 *  column defined by the row @ref _clusterbasis "clusterbasis"
 *  <tt>cb</tt> via <tt>cnext, cprev</tt>.
 *
 *  If it is necessary to update these lists, please use the functions
 *  @ref ref_row_uniform, @ref ref_col_uniform,  @ref unref_row_uniform
 *  and @ref unref_col_uniform to perform this task.
 */
struct _uniform {
  /** @brief Row @ref _clusterbasis "clusterbasis" */
  pclusterbasis rb;
  /** @brief Column @ref _clusterbasis "clusterbasis" */
  pclusterbasis cb;

  /** @brief Coupling matrix */
  amatrix S;
};

/* ------------------------------------------------------------
 * Constructors and destructors
 * ------------------------------------------------------------ */

/**
 * @brief Create a new @ref _uniform "uniform" object.
 *
 * Allocates storage for the object and sets the @ref _clusterbasis "clusterbasis"
 * pointers according to the input parameters.
 *
 * @remark Should always be matched by a call to @ref del_uniform.

 * @param rb Row @ref _clusterbasis "clusterbasis"
 * @param cb Column @ref _clusterbasis "clusterbasis"
 * @return Returns a pointer to a newly allocated @ref _uniform "uniform"
 * object.
 */
HEADER_PREFIX puniform
new_uniform(pclusterbasis rb, pclusterbasis cb);

/** @brief Deletes an @ref _uniform "uniform" object.
 *
 *  Releases the memory corresponding to the object.
 *  Reference inside the block row and block column lists are automatically
 *  removed by this function.
 *
 *  @param u Object to be deleted. */
HEADER_PREFIX void
del_uniform(puniform u);

/* ------------------------------------------------------------
 * Statistics
 * ------------------------------------------------------------ */

/** @brief Get size of a given @ref _uniform "uniform object.
 *
 *  @param u @ref _uniform "uniform" object.
 *  @returns Size of allocated storage in bytes. */
HEADER_PREFIX size_t
getsize_uniform(pcuniform u);

/* ------------------------------------------------------------
 * Simple utility functions
 * ------------------------------------------------------------ */

/**
 * @brief Sets the coupling matrix to zero.
 *
 * @param u @ref _uniform "uniform" object.
 */
HEADER_PREFIX void
clear_uniform(puniform u);

/** @brief Copy a uniform matrix @f$G@f$.
 *
 *  @param trans Set if @f$G^*@f$ should be copied instead of @f$G@f$.
 *  @param src Source matrix.
 *  @param trg target matrix. */
HEADER_PREFIX void
copy_uniform(bool trans, pcuniform src, puniform trg);

/** @brief Clone a uniform matrix.
 *
 *  @param src Source matrix.
 *  @returns Clone of the source matrix. */
HEADER_PREFIX puniform
clone_uniform(pcuniform src);

/** @brief Scale a uniform matrix.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param u Target matrix @f$S@f$, will be overwritten by @f$\alpha S@f$. */
HEADER_PREFIX void
scale_uniform(field alpha, puniform u);

/** @brief Fill a uniform matrix with random coefficients.
 *
 *  @param u Target matrix @f$S@f$. */
HEADER_PREFIX void
random_uniform(puniform u);

/* ------------------------------------------------------------
 * Access methods
 * ------------------------------------------------------------ */

/* Avoid problems with incomplete type definitions */
#if defined(CLUSTERBASIS_TYPE_COMPLETE) && !defined(UNIFORM_COMPLETE)
#define UNIFORM_COMPLETE

#ifdef __GNUC__
INLINE_PREFIX uint
getrows_uniform(pcuniform) __attribute__ ((const,unused));
INLINE_PREFIX uint
getcols_uniform(pcuniform) __attribute__ ((const,unused));
INLINE_PREFIX pamatrix
getS_uniform(puniform) __attribute__ ((unused));
#endif

/** @brief Get the number of rows of a @ref uniform matrix @f$G=V S W^*@f$.
 *
 *  @param u Matrix @f$G@f$.
 *  @return Number of rows of @f$G@f$, i.e., number of rows of @f$V@f$. */
INLINE_PREFIX uint
getrows_uniform(pcuniform u)
{
  return u->rb->t->nidx;
}

/** @brief Get the number of columns of a @ref uniform matrix @f$G=V S W^*@f$.
 *
 *  @param u Matrix @f$G@f$.
 *  @returns Number of columns of @f$G@f$, i.e., number of rows of @f$W@f$. */
INLINE_PREFIX uint
getcols_uniform(pcuniform u)
{
  return u->cb->t->nidx;
}

/** @brief Get the factor S of a @ref uniform matrix @f$G=V S W^*@f$.
 *
 * @param u Matrix @f$G@f$.
 * @returns Factor @f$S@f$. */
INLINE_PREFIX pamatrix
getS_uniform(puniform u)
{
  return &u->S;
}

#endif

#endif

/** @} */
