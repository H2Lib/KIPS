/* ------------------------------------------------------------
 * This is the file "spatialcluster.h" of the KIPS package.
 * All rights reserved, Jonas Lorenzen 2018
 * ------------------------------------------------------------ */

/** @file spatialcluster.h
 *  @author Jonas Lorenzen
 */

#ifndef SPATIALCLUSTER_H
#define SPATIALCLUSTER_H

#include "settings.h"

/** @defgroup spatialcluster spatialcluster
 *  @brief Spatial cluster tree for a given bounding box
 *  @{ */

/** @brief Data structure for a spatial cluster. */
typedef struct _spatialcluster spatialcluster;

/** @brief Pointer to a @ref spatialcluster object. */
typedef spatialcluster *pspatialcluster;

/** @brief Pointer to a constant @ref spatialcluster object. */
typedef const spatialcluster *pcspatialcluster;

/** @brief Representation of a spatial cluster.
 *
 *  A spatial cluster is an n-dimensional axis-parallel box.
 *  It may be split into two sons by bisection. */
struct _spatialcluster {
  /** @brief Spatial dimension. */
  uint dim;

  /** @brief Minimal coordinates of the bounding box. */
  preal bmin;

  /** @brief Maximal coordinates of the bounding box. */
  preal bmax;

  /** @brief Number of sons. */
  uint sons;
  
  /** @brief Pointers to sons. */
  pspatialcluster *son;
  
  /** @brief Number of clusters in entire subtree below. */
  uint desc;
  
  /** @brief Number of indices for leaf clusters. 
      Number of indices in entire subtree for non-leaf clusters.  */
  uint nidx;
  
  /** @brief Point indices for leaf clusters. */
  uint *idx;
};

/* ------------------------------------------------------------
   Constructors and destructors
   ------------------------------------------------------------ */

/** @brief Create a new @ref spatialcluster object.
 * 
 * Allocates storage for the object and sets the pointers to the sons to NULL.
 * Does not allocate storage for the index set, but sets pointer to NULL. 
 * Confer @ref init_idx_spatialcluster for initialization of the index set.
 * 
 * @remark Should always be matched by a call to @ref del_spatialcluster.
 * 
 * @param dim Spatial dimension.
 * @param bmin Minimal coordinates of the bounding box.
 * @param bmax Maximal coordinates of the bounding box.
 * @param sons Number of sons.
 * @returns New @ref spatialcluster object.*/
HEADER_PREFIX pspatialcluster
new_spatialcluster(uint dim, preal bmin, preal bmax, uint sons);

/** @brief Delete a @ref spatialcluster object.
 * 
 * Releases the storage corresponding to the @ref spatialcluster object.
 * If the cluster has sons, their storage is released too.
 * 
 * @param s @ref spatialcluster object to be deleted.*/
HEADER_PREFIX void
del_spatialcluster(pspatialcluster s);

/** @brief Update the number of indices for higher-level @ref spatialcluster objects.
 * 
 * Recursively adds alls numbers of indices of the descendants of the given 
 * @ref spatialcluster object once the index sets of the leaf clusters have been set.
 * Also updates the number of descendants.
 * 
 * @param s The @ref spatialcluster object to be updated. */
HEADER_PREFIX void
update_spatialcluster (pspatialcluster s);

/* ------------------------------------------------------------
   Bounding box
   ------------------------------------------------------------ */

/** @brief Compute the diameter of a @ref spatialcluster object in the Euclidean norm.
 * 
 * @param s Given @ref spatialcluster object.
 * @returns Euclidean norm of the diameter of @f$ s @f$.*/
HEADER_PREFIX real
diam_spatialcluster(pcspatialcluster s);

/** @brief Compute the distance of two @ref spatialcluster objects in the Euclidean norm.
 * 
 * @param s Given @ref spatialcluster object.
 * @param t Given @ref spatialcluster object.
 * @returns Euclidean norm of the distance of @f$ s @f$ and @f$ t @f$.*/
HEADER_PREFIX real
dist_spatialcluster(pcspatialcluster s, pcspatialcluster t);

/** @brief Compute the minimal Euclidean distance among all replicas of two 
 * @ref spatialcluster objects in a periodic grid.
 * 
 * The periodic grid is determined by the bounding box which contains all original 
 * @ref spatialcluster objects. This unit cell and the clusters therein are 
 * replicated in every spatial direction.
 * 
 * @param s Given @ref spatialcluster object.
 * @param t Given @ref spatialcluster object.
 * @param a Minimal coordinates of the unit cell.
 * @param b Maximal coordinates of the unit cell.*/
HEADER_PREFIX real
distPeriodic_spatialcluster (pcspatialcluster s, pcspatialcluster t, pcreal a, pcreal b);

/** @brief Compute the minimal Euclidean distance of a given point to all replicas of a 
 * @ref spatialcluster object in a periodic grid.
 * 
 * The periodic grid is determined by the bounding box which contains all original 
 * @ref spatialcluster objects. This unit cell and the clusters therein are 
 * replicated in every spatial direction. The point is required to be inside the 
 * original bounding box.
 * 
 * @param x Given point.
 * @param t Given @ref spatialcluster object.
 * @param a Minimal coordinates of the unit cell.
 * @param b Maximal coordinates of the unit cell.*/
HEADER_PREFIX real
distPeriodicPoint_spatialcluster (pcreal x, pcspatialcluster s, pcreal a, pcreal b);

/* ------------------------------------------------------------
   Statistics
   ------------------------------------------------------------ */

/** @brief Count all currently active @ref spatialcluster objects.
 * 
 * @returns Number of active @ref spatialcluster objects.*/
HEADER_PREFIX uint
getactives_spatialcluster();


/** @} */

#endif
