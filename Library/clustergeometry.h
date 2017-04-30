
/* ------------------------------------------------------------
 * This is the file "clustergeometry.h" of the KIPS package.
 * All rights reserved, Steffen Boerm 2010
 * ------------------------------------------------------------ */

#ifndef CLUSTERGEOMETRY_H
#define CLUSTERGEOMETRY_H

typedef struct _clustergeometry clustergeometry;
typedef clustergeometry *pclustergeometry;
typedef const clustergeometry *pcclustergeometry;

#include "cluster.h"
#include "settings.h"

#include <stdlib.h>

struct _clustergeometry
{
  uint dim;			/* Spatial dimension */

  uint nidx;			/* Number of indices */
  real **x;			/* Characteristic points */
  real **smin;			/* Minimal coordinates of supports */
  real **smax;			/* Maximal coordinates of supports */

  real *pmin;			/* Internal: Box for characteristic points */
  real *pmax;

  uint *perm;			/* Internal: Permutation array */
  uint *pflag;			/* Internal: Flags used for permutations */
  uint *pbuf;			/* Internal: Buffer used for permutations */

  real *buf;			/* Internal: Auxiliary storage */
};

/* ------------------------------------------------------------
   Constructors and destructors
   ------------------------------------------------------------ */

HEADER_PREFIX pclustergeometry
new_clustergeometry(uint dim, uint nidx);

HEADER_PREFIX void
del_clustergeometry(pclustergeometry cg);

/* ------------------------------------------------------------
   Building cluster trees, low-level routines
   ------------------------------------------------------------ */

/* Geometrically balanced splitting */
pcluster
splitgeometric_clustergeometry(pclustergeometry cg,
			       int levels, uint resolution,
			       uint size, uint *idx);

/* Geometrically regular splitting, requires pmin/pmax to contain
   a bounding box for the characteristic points in idx */
pcluster
splitregular_clustergeometry(pclustergeometry cg,
			     int levels, uint resolution, uint split,
			     uint size, uint *idx);

/* Cardinality-balanced splitting, requires perm to contain
   permutations of idx for all coordinate directions that put
   the coordinates into ascending order */
pcluster
splitcardinality_clustergeometry(pclustergeometry cg,
				 int levels, uint resolution,
				 uint size, uint *idx,
				 uint *perm);

/* Mixed splitting strategy: use cardinality-balanced approach for
   cardepth levels, then switch to geometric splitting */
pcluster
splitmixed_clustergeometry(pclustergeometry cg, 
			   uint cardepth, int levels, uint resolution,
			   uint size, uint *idx,
			   uint *perm);

/* ------------------------------------------------------------
   Building cluster trees, high-level routines
   ------------------------------------------------------------ */

pcluster
buildgeometric_clustergeometry(pclustergeometry cg,
			       uint levels, uint resolution,
			       uint size, uint *idx);

pcluster
buildregular_clustergeometry(pclustergeometry cg,
			     uint levels, uint resolution,
			     uint size, uint *idx);

pcluster
buildcardinality_clustergeometry(pclustergeometry cg,
				 uint levels, uint resolution,
				 uint size, uint *idx);

pcluster
buildmixed_clustergeometry(pclustergeometry cg, uint cardepth,
			   uint levels, uint resolution,
			   uint size, uint *idx);

#endif
