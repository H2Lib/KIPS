
/* ------------------------------------------------------------
 * This is the file "clusterbasis.h" of the KIPS package.
 * All rights reserved, Steffen Boerm 2009
 * ------------------------------------------------------------ */

/** @file clusterbasis.h
 *  @author Steffen B&ouml;rm
 */

#ifndef CLUSTERBASIS_H
#define CLUSTERBASIS_H

/** @defgroup clusterbasis clusterbasis
 *  @brief Representation of cluster bases for @f$\mathcal{H}^2@f$-matrices.
 *
 *  The @ref clusterbasis class represents a cluster basis
 *  @f$(V_t)_{t\in{\mathcal T}_{\mathcal{I}}}@f$, typically described
 *  by transfer matrices @f$(E_t)_{t\in\mathcal{T}_{\mathcal{I}}}@f$
 *  as @f$V_t = \sum_{t'\in\operatorname{sons}(t)} V_{t'} E_{t'}@f$
 *  for non-leaf clusters @f$t@f$.
 *
 *  @{ */

/** @brief Representation of a cluster basis. */
typedef struct _clusterbasis clusterbasis;

/** @brief Pointer to @ref clusterbasis object. */
typedef clusterbasis *pclusterbasis;

/** @brief Pointer to constant @ref clusterbasis object. */
typedef const clusterbasis *pcclusterbasis;

#include "cluster.h"
#include "amatrix.h"

/** @brief Representation of a cluster basis. */
struct _clusterbasis {
  /** @brief Corresponding cluster. */
  pccluster t;

  /** @brief Maximal rank */
  uint k;
  /** @brief Sum of ranks in entire subtree below <tt>t</tt> */
  uint ktree;

  /** @brief Leaf matrix @f$V_t@f$ */
  amatrix V;
  /** @brief Transfer matrix @f$E_t@f$ to father */
  amatrix E;

  /** @brief Number of sons, either <tt>t->sons</tt> or zero */
  uint sons;
  /** @brief Pointers to sons */
  pclusterbasis *son;
};

#define CLUSTERBASIS_TYPE_COMPLETE

/* ------------------------------------------------------------
 * Constructors and destructors
 * ------------------------------------------------------------ */

/** @brief Initialize a @ref clusterbasis object.
 *
 *  Sets up the components of the object.
 *  If <tt>t</tt> is not a leaf, the array <tt>son</tt> is allocated,
 *  otherwise it is set to null.
 *
 *  @remark Should always be matched by a call to @ref uninit_clusterbasis.
 *
 *  @param cb Object to be initialized.
 *  @param t Corresponding cluster.
 *  @returns Initialized @ref clusterbasis object. */
HEADER_PREFIX pclusterbasis
init_clusterbasis(pclusterbasis cb, pccluster t);

/** @brief Initialize a @ref clusterbasis object for a leaf.
 *
 *  Sets up the components of the object.
 *  Sets <tt>son</tt> to null.
 *  If <tt>t->sons>0</tt>, this yields a partial cluster basis.
 *
 *  @remark Should always be matched by a call to @ref uninit_clusterbasis.
 *
 *  @param cb Object to be initialized.
 *  @param t Corresponding cluster.
 *  @returns Initialized @ref clusterbasis object. */
HEADER_PREFIX pclusterbasis
init_leaf_clusterbasis(pclusterbasis cb, pccluster t);

/** @brief Uninitializes a @ref clusterbasis object.
 *
 *  Invalidates pointers, freeing corresponding storage if appropriate,
 *  and prepares the object for deletion.
 *
 *  If this @ref clusterbasis references sons, these sons
 *  are unreferences.
 *
 *  Only objects with <tt>cb->refs==0</tt> may be uninitialized.
 *
 *  @param cb Object to be uninitialized. */
HEADER_PREFIX void
uninit_clusterbasis(pclusterbasis cb);

/** @brief Create a new @ref clusterbasis object.
 *
 *  Allocates storage for the object and sets up its components.
 *
 *  @remark Should always be matched by a call to @ref del_clusterbasis.
 *  This call happens automatically if the last object referencing
 *  this object (see @ref ref_clusterbasis) is deleted.
 *
 *  @param t Corresponding cluster.
 *  @return Returns the newly created @ref clusterbasis object.
 */
HEADER_PREFIX pclusterbasis
new_clusterbasis(pccluster t);

/** @brief Create a new @ref clusterbasis object for a leaf.
 *
 *  Allocates storage for the object and sets up its components.
 *
 *  @remark Should always be matched by a call to @ref del_clusterbasis.
 *  This call happens automatically if the last object referencing
 *  this object (see @ref ref_clusterbasis) is deleted.
 *
 *  @param t Corresponding cluster.
 *  @return Returns the newly created @ref clusterbasis object.
 */
HEADER_PREFIX pclusterbasis
new_leaf_clusterbasis(pccluster t);

/** @brief Delete a @ref clusterbasis object.
 *
 *  Releases the storage corresponding to the object.
 *  If this @ref clusterbasis references sons, these sons
 *  are unreferenced.
 *
 *  Only objects with <tt>cb->refs==0</tt> may be deleted.
 *
 *  @param cb Object to be deleted. */
HEADER_PREFIX void
del_clusterbasis(pclusterbasis cb);

/* ------------------------------------------------------------
 * Low-level management
 * ------------------------------------------------------------ */

/** @brief Updates bookkeeping information, e.g., <tt>cb->ktree,</tt> for
 *  a @ref clusterbasis object after its sons have been altered.
 *
 *  @remark This function will not update the sons recursively.
 *  See @ref update_tree_clusterbasis if you really need to do this.
 *  
 *  @param cb @ref clusterbasis that will be updated. */
HEADER_PREFIX void
update_clusterbasis(pclusterbasis cb);

/** @brief Change the rank of a cluster basis and resize
 *  <tt>cb->V</tt>, while <tt>cb->son[i]->E</tt> is resized for all sons.
 *
 *  @remark <tt>cb->E</tt> is not resized, since this is in most standard
 *  algorithms a task for the father cluster.
 *
 *  @param cb Cluster basis that will be changed.
 *  @param k New rank, i.e., number of columns of <tt>V</tt> or <tt>cb->son[i]->E</tt>. */
HEADER_PREFIX void
setrank_clusterbasis(uint k, pclusterbasis cb);

/* ------------------------------------------------------------
 * Build clusterbasis based on cluster
 * ------------------------------------------------------------ */

/** @brief Construct a @ref clusterbasis from a cluster tree.
 *
 *  All ranks will be set to zero.
 *
 *  @param t Root cluster.
 *  @returns New root @ref clusterbasis object following the
 *         structure of the cluster tree. */
HEADER_PREFIX pclusterbasis
build_from_cluster_clusterbasis(pccluster t);

/* ------------------------------------------------------------
 * Statistics
 * ------------------------------------------------------------ */

/** @brief Get number of active @ref clusterbasis objects.
 *
 *  @returns Number of active @ref clusterbasis objects. */
HEADER_PREFIX uint
getactives_clusterbasis();

/** @brief Get the size of a given @ref clusterbasis object.
 *
 *  Yields the total size of this object and all its descendants.
 *
 *  @param cb @ref clusterbasis root.
 *  @returns Size of allocated storage in bytes. */
HEADER_PREFIX size_t
getsize_clusterbasis(pcclusterbasis cb);

/* ------------------------------------------------------------
 * Forward and backward transformation
 * ------------------------------------------------------------ */

/** @brief Forward transformation.
 *
 *  Compute @f$\hat x_t = V_t^* x@f$ for all elements of the cluster basis.
 *  This function also stores the permuted coefficients corresponding to
 *  the leaves of the cluster basis in the result, preparing all the
 *  necessary information for the multiplication phase realized in
 *  @ref fastaddeval_h2matrix_avector and @ref fastaddevaltrans_h2matrix_avector.
 *
 *  If <tt>cb</tt> is not a leaf, the first <tt>cb->k</tt> rows of
 *  <tt>xt</tt> are filled with @f$\hat x_t = V_t^* x@f$.
 *  The following <tt>cb->son[0]->ktree</tt> entries are filled
 *  with the coefficients for the first son, proceeding until the
 *  last <tt>cb->son[sons-1]->ktree</tt> entries are filled with the
 *  coefficients for the last son.
 *
 *  If <tt>cb</tt> is a leaf, the first <tt>cb->k</tt> rows of
 *  <tt>xt</tt> are also filled with @f$\hat x_t = V_t^* x@f$.
 *  The following <tt>cb->t->size</tt> rows are filled with the
 *  coefficients of the original vector <tt>x</tt> using the cluster
 *  numbering of <tt>cb->t</tt>.
 *
 *  @remark This function accesses the source vector <tt>x</tt>
 *  via the indices in <tt>cb->t->idx</tt>, so even if <tt>cb</tt>
 *  corresponds only to a small cluster, <tt>x</tt> usually has to
 *  be a vector corresponding to the root of the cluster tree in the
 *  original numbering.
 *  In order to work efficiently with subvectors, consider using
 *  @ref forward_nopermutation_clusterbasis_avector .
 *
 *  @param cb Cluster basis.
 *  @param x Source vector.
 *  @param xt Target vector of dimension <tt>cb->ktree</tt>, will
 *         be filled with a mix of transformed coefficients and
 *         permuted coefficients. */
HEADER_PREFIX void
forward_clusterbasis(pcclusterbasis cb, pcavector x, pavector xt);

/** @brief Backward transformation.
 *
 *  Compute @f$y \gets y + V_t \widehat y_t@f$ for all elements
 *  of the cluster basis.
 *  This function also adds permuted coefficients contained in
 *  <tt>yt</tt> to the result vector.
 *  It can be used to add the result of @ref fastaddeval_h2matrix_avector or
 *  @ref fastaddevaltrans_h2matrix_avector to the target vector.
 *
 *  The contents of <tt>yt</tt> are interpreted as in the function
 *  @ref forward_clusterbasis_avector :
 *  if <tt>cb</tt> is not a leaf, the first <tt>cb->k</tt> rows of
 *  <tt>yt</tt> are coefficients to be multiplied by @f$V_t@f$,
 *  followed by <tt>cb->son[0]->ktree</tt> coefficients for the
 *  first son, and <tt>cb->son[i]->ktree</tt> coefficients for the
 *  <tt>i</tt>-th son, until the last son with <tt>i==sons-1</tt>
 *  has been reached.
 *
 *  If <tt>cb</tt> is a leaf, the first <tt>cb->k</tt> rows of
 *  <tt>yt</tt> are also coefficients to be multiplied by @f$V_t@f$,
 *  followed by <tt>cb->t->size</tt> coefficients in cluster numbering
 *  to be added to appropriate entries of the target vector.
 *
 *  @remark This function accesses the target vector <tt>y</tt>
 *  via the indices in <tt>cb->t->idx</tt>, so even if <tt>cb</tt>
 *  corresponds only to a small cluster, <tt>y</tt> usually has to
 *  be a vector corresponding to the root of the cluster tree in
 *  the original numbering.
 *  In order to work efficiently with subvectors, consider using
 *  @ref backward_nopermutation_clusterbasis_avector .
 *
 *  @param cb Cluster basis.
 *  @param yt Source vector of dimension <tt>cb->ktree</tt>, filled
 *         with a mix of transformed coefficients and permuted coefficients.
 *         This vector will be overwritten by the function.
 *  @param y Target vector. */
HEADER_PREFIX void
backward_clusterbasis(pcclusterbasis cb, pavector yt, pavector y);

/* ------------------------------------------------------------
 * Access methods
 * ------------------------------------------------------------ */

/* Avoid problems with incomplete type definitions */
#if defined(CLUSTERBASIS_TYPE_COMPLETE) && !defined(CLUSTERBASIS_COMPLETE)
#define CLUSTERBASIS_COMPLETE

#ifdef __GNUC__
INLINE_PREFIX pamatrix
getE_clusterbasis(pclusterbasis) __attribute__ ((unused));
INLINE_PREFIX pamatrix
getV_clusterbasis(pclusterbasis) __attribute__ ((unused));
#endif

/** @brief Get the transfer matrix @f$E@f$ to the father cluster
 * of a @ref clusterbasis.
 *
 * @param cb Cluster basis.
 * @returns Transfer matrix @f$E@f$ */
INLINE_PREFIX pamatrix
getE_clusterbasis(pclusterbasis cb)
{
  return &cb->E;
}

/** @brief Get the leaf matrix @f$V@f$ of a @ref clusterbasis.
 *
 * @param cb Cluster basis.
 * @returns Leaf matrix @f$V@f$ */
INLINE_PREFIX pamatrix
getV_clusterbasis(pclusterbasis cb)
{
  return &cb->V;
}

#endif

#endif

/** @} */
