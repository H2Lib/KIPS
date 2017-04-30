
/* ------------------------------------------------------------
 * This is the file "h2matrix.h" of the KIPS package.
 * All rights reserved, Steffen Boerm 2009
 * ------------------------------------------------------------ */

/** @file h2matrix.h
 *  @author Steffen B&ouml;rm
 */

#ifndef H2MATRIX_H
#define H2MATRIX_H

/** @defgroup h2matrix h2matrix
 *  @brief Representation of an @f$\mathcal{H}^2@f$-matrix.
 *
 *  The @ref h2matrix class is used to represent @f$\mathcal{H}^2@f$-matrices
 *  with arbitrary block structures and arbitrary cluster bases.
 *  @{ */

/** @brief Representation of an @f$\mathcal{H}^2@f$-matrix. */
typedef struct _h2matrix h2matrix;

/** @brief Pointer to @ref h2matrix object. */
typedef h2matrix *ph2matrix;

/** @brief Pointer to constant @ref h2matrix object. */
typedef const h2matrix *pch2matrix;

#ifdef USE_CAIRO
#include <cairo/cairo.h>
#endif

#include "amatrix.h"
#include "block.h"
#include "uniform.h"
#include "clusterbasis.h"
#include "settings.h"

/** @brief Representation of @f$\mathcal{H}^2@f$-matrices.
 *
 *  @f$\mathcal{H}^2@f$-matrices are represented recursively:
 *  an @ref h2matrix object can be either an @ref amatrix,
 *  a @ref uniform matrix or divided into submatrices represented
 *  again by @ref h2matrix objects. */
struct _h2matrix {
  /** @brief Row cluster basis. */
  pclusterbasis rb;
  /** @brief Column cluster basis. */
  pclusterbasis cb;

  /** @brief Uniform matrix, for admissible leaves. */
  puniform u;

  /** @brief Standard matrix, for inadmissible leaves. */
  pamatrix f;

  /** @brief Submatrices. */
  ph2matrix *son;
  /** @brief Number of block rows. */
  uint rsons;
  /** @brief Number of block columns. */
  uint csons;

  /** @brief Number of descendants in matrix tree. */
  uint desc;
};

/* ------------------------------------------------------------
 * Constructors and destructors
 * ------------------------------------------------------------ */

/** @brief Create a new @ref h2matrix object.
 *
 *  Allocates storage for the object and sets the matrix pointers
 *  to <tt>NULL,</tt> representing a zero matrix.
 *
 *  @remark Should always be matched by a call to @ref del_h2matrix.
 *  This call happens automatically if the last object referencing
 *  this object (see @ref ref_h2matrix) is deleted.
 *
 *  @param rb Row cluster basis.
 *  @param cb Column cluster basis.
 *  @returns New @ref h2matrix object. */
HEADER_PREFIX ph2matrix
new_h2matrix(pclusterbasis rb, pclusterbasis cb);

/** @brief Create a new @ref h2matrix object representing a
 *  @ref uniform matrix.
 *
 *  Allocates storage for the object containing a @ref uniform
 *  object.
 *
 *  @remark Should always be matched by a call to @ref del_h2matrix.
 *  This call happens automatically if the last object referencing
 *  this object (see @ref ref_h2matrix) is deleted.
 *
 *  @param rb Row cluster basis.
 *  @param cb Column cluster basis.
 *  @returns New @ref h2matrix object containing a new @ref uniform matrix. */
HEADER_PREFIX ph2matrix
new_uniform_h2matrix(pclusterbasis rb, pclusterbasis cb);

/** @brief Create a new @ref h2matrix object representing a
 *  standard dense matrix.
 *
 *  Allocates storage for the object containing an @ref amatrix
 *  object.
 *
 *  @remark Should always be matched by a call to @ref del_h2matrix.
 *  This call happens automatically if the last object referencing
 *  this object (see @ref ref_h2matrix) is deleted.
 *
 *  @param rb Row cluster basis.
 *  @param cb Column cluster basis.
 *  @returns New @ref h2matrix object containing a new @ref amatrix. */
HEADER_PREFIX ph2matrix
new_full_h2matrix(pclusterbasis rb, pclusterbasis cb);

/** @brief Create a new @ref hmatrix object representing a
 *  subdivided matrix.
 *
 *  Allocates storage for the object representing a matrix with
 *  submatrices. The submatrices are initialized with NULL pointers
 *  and have to be set up with @ref ref_hmatrix.
 *
 *  @remark Should always be matched by a call to @ref del_h2matrix.
 *  This call happens automatically if the last object referencing
 *  this object (see @ref ref_h2matrix) is deleted.
 *
 *  @param rb Row cluster basis.
 *  @param cb Column cluster basis.
 *  @param rsons Number of block rows.
 *  @param csons Number of block columns.
 *  @returns New @ref h2matrix object with submatrices. */
HEADER_PREFIX ph2matrix
new_super_h2matrix(pclusterbasis rb, pclusterbasis cb, uint rsons, uint csons);

/** @brief Create a new @ref h2matrix object representing a zero matrix.
 *
 *  Allocates storage for the object representing a zero matrix.
 *
 *  @remark Should always be matched by a call to @ref del_h2matrix.
 *  This call happens automatically if the last object referencing
 *  this object (see @ref ref_h2matrix) is deleted.
 *
 *  @param rb Row cluster basis.
 *  @param cb Column cluster basis.
 *  @returns New @ref h2matrix object. */
HEADER_PREFIX ph2matrix
new_zero_h2matrix(pclusterbasis rb, pclusterbasis cb);

/**
 * @brief Builds a new @ref h2matrix with @ref clusterbasis <tt>rb</tt> and
 * <tt>cb</tt> and the block structure of <tt>h2</tt>.
 *
 * @param h2 Input @ref h2matrix whose structure has to be cloned.
 * @param rb Row @ref clusterbasis for the clone.
 * @param cb Column @ref clusterbasis for the clone.
 * @return New @ref h2matrix with @ref clusterbasis <tt>rb</tt> and
 * <tt>cb</tt> and the same block structure as <tt>h2</tt>.
 */
HEADER_PREFIX ph2matrix
clonestructure_h2matrix(pch2matrix h2, pclusterbasis rb, pclusterbasis cb);

/* clones a h2matrix, should be in the module h2matrix */
/**
 * @brief Builds a new @ref h2matrix with @ref clusterbasis <tt>rb</tt> and
 * <tt>cb</tt> and the block structure of <tt>h2</tt>. Copies all coupling
 * matrices @f$ S_b,\ b \in \mathcal L_{\mathcal I \times \mathcal J}^+@f$ to
 * the clone.
 *
 * @param h2 Input @ref h2matrix which has to be cloned.
 * @param rb Row @ref clusterbasis for the clone.
 * @param cb Column @ref clusterbasis for the clone.
 * @return New @ref h2matrix with @ref clusterbasis <tt>rb</tt> and
 * <tt>cb</tt> and the same data as <tt>h2</tt>.
 */
HEADER_PREFIX ph2matrix
clone_h2matrix(pch2matrix h2, pclusterbasis rb, pclusterbasis cb);

/** @brief Complete the initialisation of a @ref h2matrix object.
 *
 * Complete the initialisation of the @ref h2matrix object after all sons
 * have been initialised. 
 * The number of the descendants of the @f$\mathcal{H}^2@f$-matrix is
 * computed from the <tt>desc</tt> fields of its sons.
 * 
 * @param h2 @f$\mathcal{H}^2@f$-matrix to be completed. */
HEADER_PREFIX void
update_h2matrix(ph2matrix h2);

/** @brief Delete an @ref h2matrix object.
 *
 *  Releases the storage corresponding to the object.
 *  If this @ref h2matrix contains pointers to submatrices,
 *  the submatrices are released by @ref unref_h2matrix.
 *
 *  Only objects with <tt>h2->refs==0</tt> may be deleted.
 *
 *  @param h2 Object to be deleted. */
HEADER_PREFIX void
del_h2matrix(ph2matrix h2);

/* ------------------------------------------------------------
 * Statistics
 * ------------------------------------------------------------ */

/** @brief Get size of a given @ref h2matrix object.
 *
 *  @param h2 @f$\mathcal{H}^2@f$-matrix object.
 *  @returns Size of allocated storage in bytes. */
HEADER_PREFIX size_t
getsize_h2matrix(pch2matrix h2);

/** @brief Get total size of a given @ref h2matrix object, including
 *  cluster bases.
 *
 *  @param h2 @f$\mathcal{H}^2@f$-matrix object.
 *  @returns Size of allocated storage in bytes. */
HEADER_PREFIX size_t
gettotalsize_h2matrix(pch2matrix h2);

/** @brief Get size of the nearfield part of a given @ref h2matrix object.
 *
 *  @param h2 @f$\mathcal{H}^2@f$-matrix object.
 *  @returns Size of allocated storage for nearfield in bytes. */
HEADER_PREFIX size_t
getnearsize_h2matrix(pch2matrix h2);

/** @brief Get size of the farfield part of a given @ref h2matrix object.
 *
 *  @param h2 @f$\mathcal{H}^2@f$-matrix object.
 *  @returns Size of allocated storage for farfield in bytes. */
HEADER_PREFIX size_t
getfarsize_h2matrix(pch2matrix h2);

/* ------------------------------------------------------------
 * Simple utility functions
 * ------------------------------------------------------------ */

/** @brief Set an @ref h2matrix to zero by clearing all far- and nearfield
 *  matrices.
 *
 *  @param h2 Target matrix. */
HEADER_PREFIX void
clear_h2matrix(ph2matrix h2);

/** @brief Scale an @ref h2matrix by a factor.
 *
 * @param alpha Scaling factor @f$\alpha@f$.
 * @param h2 Target matrix @f$G@f$, will be overwritten by @f$\alpha G@f$. */
void
scale_h2matrix(field alpha, ph2matrix h2);

/** @brief Fill an @ref h2matrix with random coefficients.
 *
 * @param h2 Target matrix @f$G@f$. */
void
random_h2matrix(ph2matrix h2);

/* ------------------------------------------------------------
 * Build H^2-matrix based on block tree
 * ------------------------------------------------------------ */

/** @brief Build an @ref h2matrix object from a @ref block tree using
 *  given cluster bases.
 *
 *  @remark Submatrices for far- and nearfield leaves are created,
 *  but their coefficients are not initialized.
 *
 *  @param b Block tree.
 *  @param rb Row cluster basis.
 *  @param cb Column cluster basis.
 *  @returns New @ref h2matrix object. */
HEADER_PREFIX ph2matrix
build_from_block_h2matrix(pcblock b, pclusterbasis rb, pclusterbasis cb);

/* ------------------------------------------------------------
 * Matrix-vector multiplication
 * ------------------------------------------------------------ */

/** @brief Matrix-vector multiplication
 *  @f$y \gets y + \alpha A x@f$ or @f$y \gets y + \alpha A^* x@f$.
 *
 *  The matrix or its adjoint is multiplied by the source vector @f$x@f$,
 *  the result is scaled by @f$\alpha@f$ and added to the target vector
 *  @f$y@f$.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param h2trans Set if @f$A^*@f$ is to be used instead of @f$A@f$.
 *  @param h2 Matrix @f$A@f$.
 *  @param x Source vector @f$x@f$.
 *  @param y Target vector @f$y@f$. */
HEADER_PREFIX void
mvm_h2matrix(field alpha, bool h2trans, pch2matrix h2,
	     pcavector x, pavector y);

/** @brief Interaction phase of the matrix-vector multiplication.
 *
 *  Nearfield blocks are added directly
 *  @f$y|_{\hat t} \gets y|_{\hat t} + A|_{\hat t\times\hat s} x|_{\hat s}@f$,
 *  farfield block contributions are accumulated
 *  @f$\hat y_t \gets \hat y_t + S_{t,s} \hat x_s@f$.
 *
 *  Both <tt>xt</tt> and <tt>yt</tt> should be coefficient vectors provided by
 *  @ref new_coeffs_clusterbasis_avector, <tt>xt</tt> is usually initialized by
 *  @ref forward_clusterbasis_avector, while <tt>yt</tt> is typically added to the
 *  result using @ref backward_clusterbasis_avector.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param h2 Matrix @f$A@f$.
 *  @param xt Coefficients @f$(\hat x_s)_{s\in\mathcal{T}_{\mathcal J}}@f$
 *            of the source vector with respect to the
 *            column basis <tt>h2->cb</tt>.
 *  @param yt Coefficients @f$(\hat y_t)_{t\in\mathcal{T}_{\mathcal I}}@f$
 *            of the target vector with respect to the
 *            row basis <tt>h2->rb</tt>. */
HEADER_PREFIX void
fastaddeval_h2matrix_avector(field alpha, pch2matrix h2,
			     pcavector xt, pavector yt);

/** @brief Matrix-vector multiplication
 *  @f$y \gets y + \alpha A x@f$.
 *
 *  The matrix is multiplied by the source vector @f$x@f$, the result
 *  is scaled by @f$\alpha@f$ and added to the target vector @f$y@f$.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param h2 Matrix @f$A@f$.
 *  @param x Source vector @f$x@f$.
 *  @param y Target vector @f$y@f$. */
HEADER_PREFIX void
addeval_h2matrix(field alpha, pch2matrix h2, pcavector x, pavector y);

/** @brief Interaction phase of the adjoint matrix-vector multiplication.
 *
 *  Nearfield blocks are added directly
 *  @f$y|_{\hat s} \gets y|_{\hat s} + A|_{\hat t\times\hat s}^* x|_{\hat t}@f$,
 *  farfield block contributions are accumulated
 *  @f$\hat y_s \gets \hat y_s + S_{t,s}^* \hat x_t@f$.
 *
 *  Both <tt>xt</tt> and <tt>yt</tt> should be coefficient vectors provided by
 *  @ref new_coeffs_clusterbasis_avector, <tt>xt</tt> is usually initialized by
 *  @ref forward_clusterbasis_avector, while <tt>yt</tt> is typically added to the
 *  result using @ref backward_clusterbasis_avector.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param h2 Matrix @f$A@f$.
 *  @param xt Coefficients @f$(\hat x_t)_{t\in\mathcal{T}_{\mathcal I}}@f$
 *            of the source vector with respect to the
 *            row basis <tt>h2->rb</tt>.
 *  @param yt Coefficients @f$(\hat y_s)_{s\in\mathcal{T}_{\mathcal J}}@f$
 *            of the target vector with respect to the
 *            column basis <tt>h2->cb</tt>. */
HEADER_PREFIX void
fastaddevaltrans_h2matrix(field alpha, pch2matrix h2,
			  pcavector xt, pavector yt);

/** @brief Adjoint matrix-vector multiplication
 *  @f$y \gets y + \alpha A^* x@f$.
 *
 *  The adjoint matrix is multiplied by the source vector @f$x@f$, the
 *  result is scaled by @f$\alpha@f$ and added to the target vector @f$y@f$.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param h2 Matrix @f$A@f$.
 *  @param x Source vector @f$x@f$.
 *  @param y Target vector @f$y@f$. */
HEADER_PREFIX void
addevaltrans_h2matrix(field alpha, pch2matrix h2,
		      pcavector x, pavector y);

/** @brief Symmetric matrix-vector multiplication,
 *  @f$y \gets y + \alpha A x@f$, where @f$A@f$ is assumed to be
 *  self-adjoint and only its lower triangular part is used.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param h2 Matrix @f$A@f$.
 *  @param x Source vector @f$x@f$.
 *  @param y Target vector @f$y@f$. */
HEADER_PREFIX void
addevalsymm_h2matrix_avector(field alpha, pch2matrix h2,
			     pcavector x, pavector y);

/* ------------------------------------------------------------
 * Spectral norm
 * ------------------------------------------------------------ */

/** @brief Approximate the spectral norm @f$\|H\|_2@f$ of a matrix @f$H@f$.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$H^* H@f$ and computing the square root of
 *  the resulting eigenvalue approximation.
 *
 *  @param H2 @f$\mathcal H^2@f$ matrix @f$H@f$.
 *  @returns Approximation of @f$\|H\|_2@f$. */
HEADER_PREFIX real
norm2_h2matrix(pch2matrix H2);

/** @brief Approximate the spectral norm @f$\|A-B\|_2@f$ of the difference
 *  of two matrices @f$A@f$ and @f$B@f$.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$(A-B)^* (A-B)@f$ and computing the square root
 *  of the resulting eigenvalue approximation.
 *
 *  @param a @f$\mathcal H^2@f$ matrix @f$A@f$.
 *  @param b @f$\mathcal H^2@f$ matrix @f$B@f$.
 *  @returns Approximation of @f$\|A-B\|_2@f$. */
HEADER_PREFIX real
norm2diff_h2matrix(pch2matrix a, pch2matrix b);

/* ------------------------------------------------------------
 * Drawing
 * ------------------------------------------------------------ */

#ifdef USE_CAIRO
/**
 * @brief Draw a @ref h2matrix to a cairo surface.
 *
 * @param cr Cairo surface to be drawn to.
 * @param G The @ref h2matrix that should be drawn.
 * @param storage Flag that indicates if the storage requirements for every block
 *   should be depicted into the graphic.
 * @param levels Number of levels of the @ref h2matrix that should be drawn.
 *   If @p level == 0 holds, all levels will be drawn.
 */
HEADER_PREFIX void
draw_cairo_h2matrix(cairo_t *cr, pch2matrix G, bool storage, uint levels);
#endif

#endif

/* ------------------------------------------------------------
 * Access methods
 * ------------------------------------------------------------ */

#ifndef H2MATRIX_COMPLETE
#define H2MATRIX_COMPLETE

#ifdef __GNUC__
INLINE_PREFIX uint
getrows_h2matrix(pch2matrix) __attribute__ ((const,unused));
INLINE_PREFIX uint
getcols_h2matrix(pch2matrix) __attribute__ ((const,unused));
#endif

/** @brief Get the number of rows of an @ref h2matrix @f$G@f$.
 *
 *  @param h2 Matrix @f$G@f$.
 *  @return Number of rows of @f$G@f$. */
INLINE_PREFIX uint
getrows_h2matrix(pch2matrix h2)
{
  return h2->rb->t->size;
}

/** @brief Get the number of columns of an @ref h2matrix @f$G@f$.
 *
 *  @param h2 Matrix @f$G@f$.
 *  @returns Number of columns of @f$G@f$. */
INLINE_PREFIX uint
getcols_h2matrix(pch2matrix h2)
{
  return h2->cb->t->size;
}

#endif

/** @} */
