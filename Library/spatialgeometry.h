/* ------------------------------------------------------------
 * This is the file "spatialgeometry.h" of the KIPS package.
 * All rights reserved, Jonas Lorenzen 2018
 * ------------------------------------------------------------ */

/** @file spatialgeometry.h
 *  @author Jonas Lorenzen
 */

#ifndef SPATIALGEOMETRY_H
#define SPATIALGEOMETRY_H

#include "spatialcluster.h"

/** @defgroup spatialgeometry spatialgeometry
 *  @brief Administration tools for spatial clustering methods
 *  
 *  The @ref spatialgeometry object is used to build a spatial 
 *  cluster tree and maintaining easy access to every cluster.
 *  @{ */

/** @brief Data structure for a spatial geometry. */
typedef struct _spatialgeometry spatialgeometry;

/** @brief Pointer to a @ref spatialgeometry object. */
typedef spatialgeometry *pspatialgeometry;

/** @brief Pointer to a constant @ref spatialgeometry object. */
typedef const spatialgeometry *pcspatialgeometry;

/** @brief Representation of a spatial geometry.
 *
 * This auxiliary structure allocates storage for a list of all
 * @ref spatialcluster objects in a spatial cluster tree.
 * It also provides information about the depth and the 
 * root bounding box as well as the number of bisection steps
 * in a respective direction. */
struct _spatialgeometry {
  /** @brief Spatial dimension. */
  uint dim;

  /** @brief Minimal coordinates of the root bounding box. */
  preal bmin;

  /** @brief Maximal coordinates of the root bounding box. */
  preal bmax;

  /** @brief Number of cluster tree levels. */
  uint depth;
  
  /** @brief Numbers of bisection steps done. The first index marks 
   * the level, the second is for the coordinate number. */
  uint **splits;
  
  /** @brief List of all @ref spatialcluster objects in the cluster tree.
   * The first index marks the level, the second counts the number of 
   * the cluster in lexicographical order. */
  pspatialcluster **s;
};

/* ------------------------------------------------------------
   Constructors and destructors
   ------------------------------------------------------------ */

/** @brief Create a new empty @ref spatialgeometry object.
 * 
 * Only allocates storage for the object itself and sets up the
 * bounding box. For a complete initialization use @ref init_spatialgeometry.
 * 
 * @remark Should always be matched by a call to @ref del_spatialgeometry.
 * 
 * @param dim Spatial dimension.
 * @param bmin Minimal coordinates of the bounding box.
 * @param bmax Maximal coordinates of the bounding box.
 * @returns New @ref spatialgeometry object.*/
HEADER_PREFIX pspatialgeometry
new_spatialgeometry(uint dim, preal bmin, preal bmax);

/** @brief Delete a @ref spatialgeometry object.
 * 
 * Releases the storage corresponding to the @ref spatialgeometry object.
 * 
 * @param sg @ref spatialgeometry object to be deleted.*/
HEADER_PREFIX void
del_spatialgeometry(pspatialgeometry sg);

/* ------------------------------------------------------------
   Adaptive spatial cluster tree
   ------------------------------------------------------------ */

/** @brief Build a @ref spatialcluster tree using a @ref spatialgeometry
 * object and adaptive clustering.
 * 
 * Builds an adaptive @ref spatialcluster tree from the bounding box
 * given in a @ref spatialgeometry object. Bisects every @ref spatialcluster
 * in the coordinate of maximal extension, until either the given maximal
 * depth is reached or every cluster on the current level has a diameter 
 * not greater than the given maximal diameter. In particular, all 
 * @ref spatialcluster objects on one level have congruent bounding boxes.
 * Information about the bisection steps and the pointers to the 
 * @ref spatialcluster objects are stored in the respective 
 * @ref spatialgeometry object.
 * 
 * @param maxdepth Maximal depth of the resulting @ref spatialcluster tree.
 * @param maxdiam Maximal diameter for a leaf @ref spatialcluster. Note 
 * that the algorith will stop when either maxdepth OR maxdiam ar reached. 
 * @param sg The @ref spatialgeometry object containing the root bounding 
 * box. Will be updated with bookkeeping information.
 * @returns The root @ref spatialcluster.*/
HEADER_PREFIX pspatialcluster
init_spatialgeometry(uint maxdepth, real maxdiam, pspatialgeometry sg);

/** @brief Find a @ref spatialcluster containing a certain point inside a
 * given @ref spatialgeometry.
 * 
 * Finds the respective @ref spatialcluster object with bounding box 
 * containing the given point in three-dimensional space on a given 
 * level of a @ref spatialcluster tree administrated by a 
 * @ref spatialgeometry object.
 * 
 * @param l Level of the @ref spatialcluster tree on which the 
 * corresponding cluster shall be found.
 * @param x Point to be located in the @ref spatialcluster tree.
 * @param sg The @ref spatialgeometry to look in.
 * @param nr The lexicographical number of the desired @ref spatialcluster 
 * object on the given level. (To be set by calling this function.)
 * @returns The desired @ref spatialcluster object.*/
HEADER_PREFIX pspatialcluster
findCluster_spatialgeometry(uint l, pcreal x, pcspatialgeometry sg, uint *nr);

/** @brief Distribute the indices of a vector array to the corresponding 
 * leaf @ref spatialcluster objects of a @ref spatialgeometry.
 * 
 * Allocates storage for the index sets of all leaf @ref spatialcluster 
 * objects in the given @ref spatialgeometry and initializes them with 
 * the indices of those points that are contained in the respective 
 * @ref spatialcluster.
 * 
 * @param nidx Number of vectors to be distributed.
 * @param x Vector array to be distributed.
 * @param sg The @ref spatialgeometry object to distribute in. */
HEADER_PREFIX void
initPoints_spatialgeometry (uint nidx, pcreal *x, pspatialgeometry sg);

/* ------------------------------------------------------------
   Derived quantities
   ------------------------------------------------------------ */

/** @brief Determine the total number of @ref spatialcluster objects on 
 * one level of a given @ref spatialgeometry.
 * 
 * @param l Level of the @ref spatialcluster tree.
 * @param sg The @ref spatialgeometry to look in.
 * @returns Number of @ref spatialcluster objects on that level. */
HEADER_PREFIX uint
countonlevel_spatialgeometry (uint l, pcspatialgeometry sg);

/** @brief Determine the number of @ref spatialcluster objects in one 
 * spatial direction on one level of a given @ref spatialgeometry.
 * 
 * @param l Level of the @ref spatialcluster tree.
 * @param d Spatial direction.
 * @param sg The @ref spatialgeometry object to look in.
 * @returns Number of @ref spatialcluster objects in the given direction
 * on that level. */
HEADER_PREFIX uint
countindirection_spatialgeometry (uint l, uint d, pcspatialgeometry sg);


/** @} */

#endif
