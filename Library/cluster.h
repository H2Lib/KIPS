
/* ------------------------------------------------------------
 * This is the file "cluster.h" of the KIPS package.
 * All rights reserved, Steffen Boerm 2010
 * ------------------------------------------------------------ */

#ifndef CLUSTER_H
#define CLUSTER_H

typedef struct _cluster cluster;
typedef cluster *pcluster;
typedef const cluster *pccluster;

#include "settings.h"
#include <stdlib.h>

struct _cluster
{
  uint size;			/* Number of indices */
  uint *idx;			/* Index array */

  uint sons;			/* Number of sons */
  pcluster *son;		/* Pointers to sons */

  uint dim;			/* Dimension of bounding box */
  real *bmin;			/* Minimal coordinates */
  real *bmax;			/* Maximal coordinates */

  uint desc;			/* Number of descendants */
};

/* ------------------------------------------------------------
   Constructors and destructors
   ------------------------------------------------------------ */

HEADER_PREFIX pcluster
new_cluster(uint size, uint *idx, uint sons, uint dim);

HEADER_PREFIX void
update_cluster(pcluster t);

HEADER_PREFIX void
del_cluster(pcluster t);

HEADER_PREFIX void
setsons_cluster(pcluster t, uint sons);

/* ------------------------------------------------------------
   Bounding box
   ------------------------------------------------------------ */

HEADER_PREFIX real
diam_cluster(pccluster t);

HEADER_PREFIX real
dist_cluster(pccluster t, pccluster s);

/* ------------------------------------------------------------
   Statistics
   ------------------------------------------------------------ */

HEADER_PREFIX uint
getactives_cluster();

HEADER_PREFIX size_t
getsize_cluster(pccluster t);

HEADER_PREFIX uint
getdepth_cluster(pccluster t);

#endif
