
/* ------------------------------------------------------------
 * This is the file "kernelmatrix.h" of the KIPS package.
 * All rights reserved, Steffen Boerm 2018
 * ------------------------------------------------------------ */

/** @file kernelmatrix.h
 *  @author Steffen B&ouml;rm
 */

#ifndef KERNELMATRIX_H
#define KERNELMATRIX_H

#include "settings.h"
#include "h2matrix.h"
#include "spatialcluster.h"

/** @defgroup kernelmatrix kernelmatrix
 *  @brief Approximation of kernel matrices.
 *  @{ */

/** @brief Data required to approximate a kernel matrix. */
typedef struct _kernelmatrix kernelmatrix;

/** @brief Pointer to a @ref kernelmatrix object. */
typedef kernelmatrix *pkernelmatrix;

/** @brief Pointer to a constant @ref kernelmatrix object. */
typedef const kernelmatrix *pckernelmatrix;

/** @brief Pointer to a kernel function. */
typedef real (*kernel) (pcreal xx, pcreal yy, uint i, uint j, void *data);

/** @brief Pointer to the gradient of a kernel function. */
typedef void (*gradient) (pcreal xx, pcreal yy, uint i, uint j, void *data, 
                          pfield f);

/** @brief Representation of a kernel matrix an its approximation.
 *
 *  A kernel matrix is a matrix with entries of the form
 *  @f$g_{ij} = k(x_i,x_j)@f$, where @f$k@f$ is a kernel functions
 *  and @f$(x_i)_{i\in\Idx}@f$ are points in a suitable space.
 *
 *  If the kernel functions is locally smooth, it can be approximated
 *  by interpolation, and this gives rise to blockwise low-rank
 *  approximations. */
struct _kernelmatrix {
  /** @brief Spatial dimension. */
  uint dim;

  /** @brief Kernel function. */
  kernel g;
  
  /** @brief Gradient of the kernel function. */
  gradient f;

  /** @brief Data for the kernel function. */
  void *data;

  /** @brief Number of points. */
  uint points;

  /** @brief Coordinates of points. */
  real **x;
  
  /** @brief Type flags for each point, e. g. molecule number. 
   *        Zero is used for interpolation points exclusively. */
  uint *type;

  /** @brief Interpolation order (i.e., number of interpolation points). */
  uint m;

  /** @brief Interpolation points for the reference interval @f$[-1,1@f$. */
  real *xi_ref;
};

/** @brief Create an empty @ref kernelmatrix object.
 *
 *  @remark Does not allocate storage for point coordinates and type flags.
 *
 *  @param dim Spatial dimension.
 *  @param m Interpolation order.
 *  @returns New object. */
HEADER_PREFIX pkernelmatrix
new_kernelmatrix(uint dim, uint m);

/** @brief Delete a @ref kernelmatrix object.
 *
 *  @remark Does not free storage for point coordinates and type flags.
 *
 *  @param km Object to be deleted. */
HEADER_PREFIX void
del_kernelmatrix(pkernelmatrix km);

/** @brief Update a @ref kernelmatrix object with new points.
 *
 *  @remark Does not free the storage of previous coordinate or flag arrays.
 *
 *  @param points New number of points.
 *  @param x Location of new points.
 *  @param type Type flags of new points.
 *  @param km Object to be updated. */
HEADER_PREFIX void
update_kernelmatrix (uint points, preal *x, uint *type, pkernelmatrix km);

HEADER_PREFIX void
fillN_kernelmatrix (bool gradient, const uint *ridx, const uint *cidx, 
                    pckernelmatrix km, pamatrix N);

HEADER_PREFIX void
fillS_kernelmatrix(pcspatialcluster rc, pcspatialcluster cc, pckernelmatrix km, 
                   pamatrix S);

HEADER_PREFIX void
fillV_kernelmatrix (bool gradient, pcspatialcluster tc, pckernelmatrix km, 
                    pamatrix V);

HEADER_PREFIX void
fillE_kernelmatrix (pcspatialcluster sc, pcspatialcluster fc, pckernelmatrix km, 
                   pamatrix E);

/** @brief Fill a @ref clusterbasis using interpolation.
 *
 *  @param gradient Type flag for a gradient cluster basis.
 *  @param km Description of the kernel matrix.
 *  @param cb Cluster basis to be filled. */
HEADER_PREFIX void
fill_clusterbasis_kernelmatrix(bool gradient, pckernelmatrix km, pclusterbasis cb);

/** @brief Fill a @ref h2matrix using interpolation.
 *
 *  @param gradient Type flag for a gradient kernel matrix.
 *  @param km Description of the kernel matrix.
 *  @param G Matrix to be filled. */
HEADER_PREFIX void
fill_h2matrix_kernelmatrix(bool gradient, pckernelmatrix km, ph2matrix G);

/** @brief Update a @ref h2matrix object based on a @ref kernelmatrix object.
 *
 *  The algorithm will preserve the overall structure of the h2matrix while 
 *  changing leaf matrices and nearfield matrices based on new points given 
 *  by the kernelmatrix.
 *
 *  @param gradient Type flag for a gradient kernel matrix.
 *  @param km The @ref kernelmatrix object containing the new points.
 *  @param G Matrix to be updated. */
HEADER_PREFIX void
update_h2matrix_kernelmatrix(bool gradient, pckernelmatrix km, ph2matrix G);

/** @} */

#endif
