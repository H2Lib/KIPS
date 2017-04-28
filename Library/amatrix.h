
/* ------------------------------------------------------------
 * This is the file "amatrix.h" of the H2Lib package.
 * All rights reserved, Steffen Boerm 2009
 * ------------------------------------------------------------ */

/** @file amatrix.h
 *  @author Steffen B&ouml;rm
 */

#ifndef AMATRIX_H
#define AMATRIX_H

/** @defgroup amatrix amatrix
 *  @brief Representation of a matrix as an array in column-major order.
 *
 *  The @ref amatrix class is used to handle standard linear algebra
 *  operations like matrix multiplication, factorization, solving
 *  linear systems or eigenvalue problems.
 *  @{ */

/** @brief Representation of a matrix as a column-order array. */
typedef struct _amatrix amatrix;

/** @brief Pointer to @ref amatrix object. */
typedef amatrix *pamatrix;

/** @brief Pointer to constant @ref amatrix object. */
typedef const amatrix *pcamatrix;

#include <assert.h>
#include <math.h>
#include <stdlib.h>

#include "basic.h"
#include "settings.h"
#include "avector.h"
#include "realavector.h"

/** @brief Representation of a matrix as an array in column-major order. */
struct _amatrix {
  /** @brief Matrix coefficients in column-major order, i.e., @f$a_{ij}@f$ corresponds to `a[i+j*ld]`.  */
  field *a;

  /** @brief Leading dimension, i.e., increment used to switch from one column to the next.  */
  uint ld;

  /** @brief Number of rows. */
  uint rows;
  /** @brief Number of columns.  */
  uint cols;

  /** @brief Points to owner of coefficient storage if this is a submatrix. */
  void *owner;
};

/* ------------------------------------------------------------ *
 * Constructors and destructors                                 *
 * ------------------------------------------------------------ */

/** @brief Initialize an @ref amatrix object.
 *
 *  Sets up the components of the object and allocates storage for
 *  the coefficient array.
 *
 *  @remark Should always be matched by a call to @ref uninit_amatrix.
 *
 *  @param A Object to be initialized.
 *  @param rows Number of rows.
 *  @param cols Number of columns.
 *  @returns Initialized @ref amatrix object. */
HEADER_PREFIX pamatrix
init_amatrix(pamatrix A, uint rows, uint cols);

/** @brief Initialize an @ref amatrix object to represent a submatrix.
 *
 *  Sets up the components of the object and uses part of the storage
 *  of another @ref amatrix for the coefficient array, leading to a new
 *  matrix representing a submatrix of the source.
 *  Changes to the coefficients of the new matrix also change
 *  coefficients of the source matrix.
 *
 *  @remark Should always be matched by a call to @ref uninit_amatrix that
 *  will <em>not</em> release the coefficient storage.
 *
 *  @param A Object to be initialized.
 *  @param src Source matrix.
 *  @param rows Number of rows.
 *  @param roff Row offset, should satisfy <tt>rows+roff<=src->rows</tt>.
 *  @param cols Number of columns.
 *  @param coff Column offset, should satisfy <tt>cols+coff<=src->cols</tt>.
 *  @returns Initialized @ref amatrix object. */
HEADER_PREFIX pamatrix
init_sub_amatrix(pamatrix A, pamatrix src, uint rows, uint roff, uint cols,
    uint coff);

/** @brief Uninitialize an @ref amatrix object.
 *
 *  Invalidates pointers, freeing corresponding storage if appropriate,
 *  and prepares the object for deletion.
 *
 *  @param A Object to be uninitialized. */
HEADER_PREFIX void
uninit_amatrix(pamatrix A);

/** @brief Create a new @ref amatrix object.
 *
 *  Allocates storage for the object and sets up its components.
 *
 *  @remark Should always be matched by a call to @ref del_amatrix.
 *
 *  @param rows Number of rows.
 *  @param cols Number of columns.
 *  @returns New @ref amatrix object. */
HEADER_PREFIX pamatrix
new_amatrix(uint rows, uint cols);

/** @brief Create a new @ref amatrix object representing a submatrix.
 *
 *  Allocates storage for the object, but uses part of the storage
 *  of another @ref amatrix object to keep the coefficients.
 *  Since the leading dimension of the source matrix is used, this
 *  allows us to work with a rectangular submatrix of the original
 *  matrix.
 *
 *  @remark Should always be matched by a call to @ref del_amatrix that
 *  will <em>not</em> release the coefficient storage.
 *
 *  @param src Source matrix.
 *  @param rows Number of rows.
 *  @param roff Row offset, should satisfy <tt>rows+roff<=src->rows</tt>.
 *  @param cols Number of columns.
 *  @param coff Column offset, should satisfy <tt>cols+coff<=src->cols</tt>.
 *  @returns New @ref amatrix object. */
HEADER_PREFIX pamatrix
new_sub_amatrix(pamatrix src, uint rows, uint roff, uint cols, uint coff);

/** @brief Delete an @ref amatrix object.
 *
 *  Releases the storage corresponding to the object.
 *
 *  @attention Make sure that there are no submatrix objects referring
 *  to this object, since they will otherwise keep using pointers
 *  to invalid storage.
 *
 *  @param A Object to be deleted. */
HEADER_PREFIX void
del_amatrix(pamatrix A);

/** @brief Change the dimensions of an @ref amatrix object without
 *  preserving its coefficients.
 *
 *  Allocates new storage for the coefficients and releases the
 *  old storage.
 *
 *  The matrix entries are not preserved.
 *
 *  @attention Make sure that there are no submatrix objects referring
 *  to this object, since they might otherwise keep using pointers
 *  to invalid storage.
 *
 *  @param A Matrix to be resized.
 *  @param rows New number of rows.
 *  @param cols New number of columns. */
HEADER_PREFIX void
resize_amatrix(pamatrix A, uint rows, uint cols);

/** @brief Change the dimensions of an @ref amatrix object while
 *  preserving as many of its coefficients as possible.
 *
 *  Allocates new storage for the coefficients and copies as many
 *  of the old coefficients into it.
 *  If there are more rows or columns in the new matrix, the additional
 *  coefficients are left unintialized.
 *
 *  @attention Make sure that there are no submatrix objects referring
 *  to this object, since they might otherwise keep using pointers
 *  to invalid storage.
 *
 *  @param A Matrix to be resized.
 *  @param rows New number of rows.
 *  @param cols New number of columns. */
HEADER_PREFIX void
resizecopy_amatrix(pamatrix A, uint rows, uint cols);

/* ------------------------------------------------------------
 * Access methods
 * ------------------------------------------------------------ */

#ifdef __GNUC__
INLINE_PREFIX field
getentry_amatrix(pcamatrix, uint, uint) __attribute__ ((const,unused));
INLINE_PREFIX void
setentry_amatrix(pamatrix, uint, uint, field) __attribute__((unused));
INLINE_PREFIX field
addentry_amatrix(pamatrix, uint, uint, field) __attribute__((unused));
#endif

/** @brief Read a matrix entry @f$a_{ij}@f$.
 *
 *  @param A Matrix @f$A@f$.
 *  @param row Row index @f$i@f$.
 *  @param col Column index @f$j@f$.
 *  @returns Matrix entry @f$a_{ij}@f$. */
INLINE_PREFIX field
getentry_amatrix(pcamatrix A, uint row, uint col)
{
  longindex ldA = A->ld;
#ifdef FULL_DEBUG
  assert(row < A->rows);
  assert(col < A->cols);
#endif

  return A->a[row + col * ldA];
}

/** @brief Set a matrix entry, @f$a_{ij}\gets x@f$.
 *
 *  @param A Matrix @f$A@f$.
 *  @param row Row index @f$i@f$.
 *  @param col Column index @f$j@f$.
 *  @param x New value of @f$a_{ij}@f$. */
INLINE_PREFIX void
setentry_amatrix(pamatrix A, uint row, uint col, field x)
{
  longindex ldA = A->ld;
#ifdef FULL_DEBUG
  assert(row < A->rows);
  assert(col < A->cols);
#endif

  A->a[row + col * ldA] = x;
}

/** @brief Add to a matrix entry, @f$a_{ij} \gets a_{ij} + x@f$.
 *
 *  @param A Matrix @f$A@f$.
 *  @param row Row index @f$i@f$.
 *  @param col Column index @f$j@f$.
 *  @param x Summand.
 *  @returns New value of @f$a_{ij}@f$. */
INLINE_PREFIX field
addentry_amatrix(pamatrix A, uint row, uint col, field x)
{
  longindex ldA = A->ld;
#ifdef FULL_DEBUG
  assert(row < A->rows);
  assert(col < A->cols);
#endif

  return (A->a[row + col * ldA] += x);
}

/* ------------------------------------------------------------
 * Statistics
 * ------------------------------------------------------------ */

/** @brief Get number of currently initialized @ref amatrix objects.
 *
 *  Calls to initialization functions like @ref init_amatrix and
 *  constructors like @ref new_amatrix increase an internal counter,
 *  while @ref uninit_amatrix and @ref del_amatrix decrease it.
 *
 *  @remark Use this function to check whether a program correctly cleans
 *  up temporary variables.
 *
 *  @returns Number of currently initialized @ref amatrix objects. */
HEADER_PREFIX uint
getactives_amatrix();

/** @brief Get size of a given @ref amatrix object.
 *
 *  Computes the size of the @ref amatrix object and the storage
 *  allocated for the coefficients.
 *  If the object uses the coefficients of another object
 *  (e.g., if it was created using @ref new_sub_amatrix), no coefficient
 *  storage is added.
 *
 *  @param a Matrix object.
 *  @returns Size of allocated storage in bytes. */
HEADER_PREFIX size_t
getsize_amatrix(pcamatrix a);

/** @brief Get heap size of a given @ref amatrix object.
 *
 *  Computes the size of storage allocated for the coefficients
 *  on the heap, but not for the @ref amatrix object itself.
 *  If the object uses the coefficients of another object
 *  (e.g., if it was created using @ref new_sub_amatrix), no storage
 *  is required.
 *
 *  @param a Matrix object.
 *  @returns Size of allocated heap storage in bytes. */
HEADER_PREFIX size_t
getsize_heap_amatrix(pcamatrix a);

/* ------------------------------------------------------------
 * Simple utility functions
 * ------------------------------------------------------------ */

/** @brief Set a matrix to zero.
 *
 *  @param A Target matrix. */
HEADER_PREFIX void
clear_amatrix(pamatrix A);

/** @brief Set the lower triangular part of a matrix to zero.
 *
 *  @param A Target matrix.
 *  @param strict If set, only the <em>strict</em> lower triangular
 *     part should be cleared. */
HEADER_PREFIX void
clear_lower_amatrix(bool strict, pamatrix A);

/** @brief Set the upper triangular part of a matrix to zero.
 *
 *  @param A Target matrix.
 *  @param strict If set, only the <em>strict</em> upper triangular
 *     part should be cleared. */
HEADER_PREFIX void
clear_upper_amatrix(bool strict, pamatrix A);

/** @brief Set a matrix to identity.
 *
 *  Sets all diagonal entries to one and all off-diagonal entries to zero.
 *  For square matrices, this yields an identity matrix.
 *
 *  @param A Target matrix. */
HEADER_PREFIX void
identity_amatrix(pamatrix A);

/** @brief Fill a matrix with random values.
 *
 *  @param A Target matrix. */
HEADER_PREFIX void
random_amatrix(pamatrix A);

/** @brief Fill a square matrix with random values and ensure that it is invertible
 *
 *  First the matrix is filled with random values, then the diagonal
 *  elements are set to @f$a_{ii} \gets \alpha \sum_{j=1}^n |a_{ij}|@f$,
 *  ensuring that @f$A@f$ is diagonal-dominant and therefore invertible
 *  if @f$\alpha > 1@f$.
 *
 *  @param A Target matrix.
 *  @param alpha Diagonal factor @f$\alpha@f$, should be greater than one. */
HEADER_PREFIX void
random_invertible_amatrix(real alpha, pamatrix A);

/** @brief Fill a matrix with random values and ensure that it is self-adjoint.
 *
 *  @param A Target matrix. */
HEADER_PREFIX void
random_selfadjoint_amatrix(pamatrix A);

/** @brief Fill a matrix with random values and ensure that it is positive definite.
 *
 *  First the matrix is filled with random values, ensuring that it becomes
 *  self-adjoint.
 *  Then the diagonal elements are set to
 *  @f$a_{ii} \gets \alpha \sum_{j=1}^n |a_{ij}|@f$, ensuring that
 *  @f$A@f$ is diagonal-dominant and therefore positiv definite
 *  if @f$\alpha > 1@f$.
 *
 *  @param A Target matrix.
 *  @param alpha Diagonal factor @f$\alpha@f$, should be greater than one. */
HEADER_PREFIX void
random_spd_amatrix(real alpha, pamatrix a);

/** @brief Copy a matrix into another matrix, @f$B \gets A@f$ or
 *  @f$B \gets A^*@f$.
 *
 *  The numbers of rows and columns have to match.
 *
 *  @param transA Set if @f$A^*@f$ is to be used instead of @f$A@f$.
 *  @param A Source matrix.
 *  @param B Target matrix. */
HEADER_PREFIX void
copy_amatrix(bool transA, pcamatrix A, pamatrix B);

/** @brief Create a duplicate of an existing @ref amatrix.
 *
 *  @param src Matrix to be duplicated.
 *  @returns Copy of <tt>src</tt>. */
HEADER_PREFIX pamatrix
clone_amatrix(pcamatrix src);

/** @brief Copy a matrix into another matrix, @f$B \gets A@f$ or
 *  @f$B \gets A^*@f$.
 *
 *  If @f$B@f$ is smaller than @f$A@f$ (or @f$A^*@f$), only the upper
 *  left part of @f$A@f$ is copied.
 *  If @f$A@f$ (or @f$A^*@f$) is smaller than @f$B@f$, only the upper
 *  left part of @f$B@f$ is filled.
 *
 *  @param transA Set if @f$A^*@f$ is to be used instead of @f$A@f$.
 *  @param A Source matrix.
 *  @param B Target matrix. */
HEADER_PREFIX void
copy_sub_amatrix(bool transA, pcamatrix A, pamatrix B);

/** @brief Print a matrix.
 *
 *  @param A Matrix object. */
HEADER_PREFIX void
print_amatrix(pcamatrix A);

/** @brief Print a matrix in Matlab format.
 *
 *  @param A Matrix object. */
HEADER_PREFIX void
print_matlab_amatrix(pcamatrix A);

/** @brief Check whether a matrix @f$A@f$ or its adjoint @f$A^*@f$ is isometric.
 *
 *  Compute either @f$I-A^*A@f$ or @f$I-AA^*@f$.
 *  If the return value is small, @f$A@f$ or @f$A^*@f$ are isometric
 *  matrices.
 *  For square matrices, this also means that they are orthogonal.
 *
 *  @param transA Set if @f$A^*@f$ is to be checked, otherwise @f$A@f$ is used.
 *  @param A Matrix @f$A@f$.
 *  @returns @f$\|I-A^*A\|_F@f$ if <tt>transA==false</tt> and @f$\|I-AA^*\|_F@f$ otherwise. */
HEADER_PREFIX real
check_ortho_amatrix(bool transA, pcamatrix A);

/* ------------------------------------------------------------
 * Basic linear algebra
 * ------------------------------------------------------------ */

/** @brief Scale a matrix @f$A@f$ by a factor @f$\alpha@f$,
 *  @f$A \gets \alpha A@f$.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param A Target matrix @f$A@f$. */
HEADER_PREFIX void
scale_amatrix(field alpha, pamatrix A);

/** @brief Compute the Frobenius inner product
 *  @f$\langle A, B \rangle_F@f$ of two matrices @f$A@f$ and @f$B@f$.
 *
 *  The Frobenius inner product is given by
 *  @f$\langle A, B \rangle_F = \sum_{i,j} \bar a_{ij} b_{ij}@f$.
 *
 *  @param A Matrix @f$A@f$.
 *  @param B Matrix @f$B@f$.
 *  @returns Frobenius inner product @f$\langle A, B \rangle_F@f$. */
HEADER_PREFIX field
dotprod_amatrix(pcamatrix A, pcamatrix B);

/** @brief Approximate the spectral norm @f$\|A\|_2@f$ of a matrix @f$A@f$.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$A^* A@f$ and computing the square root of
 *  the resulting eigenvalue approximation.
 *
 *  @param A Dense matrix @f$A@f$.
 *  @returns Approximation of @f$\|A\|_2@f$. */
HEADER_PREFIX real
norm2_amatrix(pcamatrix A);

/** @brief Compute the Frobenius norm @f$\|A\|_F@f$ of a matrix @f$A@f$.
 *
 *  The Frobenius norm is given by
 *  @f$\|A\|_F = \left( \sum_{i,j} |a_{ij}|^2 \right)^{1/2}@f$.
 *
 *  @param A Matrix @f$A@f$.
 *  @returns Frobenius norm @f$\|A\|_F@f$. */
HEADER_PREFIX real
normfrob_amatrix(pcamatrix A);

/** @brief Approximate the spectral norm @f$\|A-B\|_2@f$ of the difference
 *  of two matrices @f$A@f$ and @f$B@f$.
 *
 *  The spectral norm is approximated by applying a few steps of the power
 *  iteration to the matrix @f$(A-B)^* (A-B)@f$ and computing the square root
 *  of the resulting eigenvalue approximation.
 *
 *  @param A Dense matrix @f$A@f$.
 *  @param B Dense matrix @f$B@f$.
 *  @returns Approximation of @f$\|A-B\|_2@f$. */
HEADER_PREFIX real
norm2diff_amatrix(pcamatrix A, pcamatrix B);

/** @brief Multiply a matrix @f$A@f$ by a vector @f$x@f$,
 *  @f$y \gets y + \alpha A x@f$.
 *
 *  The matrix @f$A@f$ is multiplied by the source vector @f$x@f$,
 *  the result is scaled by @f$\alpha@f$ and added to the
 *  target vector @f$y@f$.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param A Matrix @f$A@f$.
 *  @param x Source vector @f$x@f$.
 *  @param y Target vector @f$y@f$. */
HEADER_PREFIX void
addeval_amatrix(field alpha, pcamatrix A,
		pcavector x, pavector y);

/** @brief Multiply the adjoint of a matrix @f$A@f$ by a vector @f$x@f$,
 *  @f$y \gets y + \alpha A^* x@f$.
 *
 *  The adjoint @f$A^*@f$ is multiplied by the source vector @f$x@f$,
 *  the result is scaled by @f$\alpha@f$ and added to the
 *  target vector @f$y@f$.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param A Matrix @f$A@f$.
 *  @param x Source vector @f$x@f$.
 *  @param y Target vector @f$y@f$. */
HEADER_PREFIX void
addevaltrans_amatrix_avector(field alpha, pcamatrix A,
			     pcavector x, pavector y);

/** @brief Multiply a matrix @f$A@f$ or its adjoint @f$A^*@f$ by a
 *  vector, @f$y \gets y + \alpha A x@f$ or @f$y \gets y + \alpha A^* x@f$.
 *
 *  The matrix or its adjoint is multiplied by the source vector @f$x@f$,
 *  the result is scaled by @f$\alpha@f$ and added to the target vector
 *  @f$y@f$.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param transA Set if @f$A^*@f$ is to be used instead of @f$A@f$.
 *  @param A Matrix @f$A@f$.
 *  @param x Source vector @f$x@f$.
 *  @param y Target vector @f$y@f$. */
HEADER_PREFIX void
mvm_amatrix_avector(field alpha, bool transA, pcamatrix a,
		    pcavector x, pavector y);

/** @brief Add two matrices,
 *  @f$B \gets B + \alpha A@f$ or @f$B \gets B + \alpha A^*@f$.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param transA Set if @f$A^*@f$ is to be added instead of @f$A@f$.
 *  @param A Source matrix @f$A@f$.
 *  @param B Target matrix @f$B@f$. */
HEADER_PREFIX void
add_amatrix(field alpha, bool transA, pcamatrix A, pamatrix B);

/** @brief Multiply two matrices,
 *  @f$C \gets C + \alpha A B@f$, @f$C \gets C + \alpha A^* B@f$,
 *  @f$C \gets C + \alpha A B^*@f$ or @f$C \gets C + \alpha A^* B^*@f$.
 *
 *  The matrices @f$A@f$ (or @f$A^*@f$) and @f$B@f$ (or @f$B^*@f$)
 *  are multiplied, the result is scaled by @f$\alpha@f$ and added
 *  to the target matrix @f$C@f$.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param transA Set if @f$A^*@f$ is to be used instead of @f$A@f$.
 *  @param A Left factor @f$A@f$.
 *  @param transB Set if @f$B^*@f$ is to be used instead of @f$B@f$.
 *  @param B Right factor @f$B@f$.
 *  @param C Target matrix @f$C@f$. */
HEADER_PREFIX void
addmul_amatrix(field alpha, bool transA, pcamatrix A,
	       bool transB, pcamatrix B, pamatrix C);

/** @} */

#endif
