
/* ------------------------------------------------------------
   This is the file "cluster.h" of the H2Lib package.
   All rights reserved, Steffen Boerm 2010
   ------------------------------------------------------------ */

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

/* ------------------------------------------------------------
   Hierarchical iterator
   ------------------------------------------------------------ */

HEADER_PREFIX void
iterate_cluster(pccluster t, uint tname,
		void (*pre)(pccluster t, uint tname, void *data),
		void (*post)(pccluster t, uint tname, void *data),
		void *data);

HEADER_PREFIX void
iterate_parallel_cluster(pccluster t, uint tname,
			 uint pardepth,
			 void (*pre)(pccluster t, uint tname, void *data),
			 void (*post)(pccluster t, uint tname, void *data),
			 void *data);

/* ------------------------------------------------------------
   Enumeration
   ------------------------------------------------------------ */

HEADER_PREFIX pcluster *
enumerate_cluster(pcluster t);

/* ------------------------------------------------------------
 * File I/O
 * ------------------------------------------------------------ */

#ifdef USE_NETCDF
HEADER_PREFIX void
write_cdf_cluster(pccluster t, const char *name);

HEADER_PREFIX void
write_cdfpart_cluster(pccluster t, int nc_file, const char *prefix);

HEADER_PREFIX pcluster
read_cdf_cluster(const char *name);

HEADER_PREFIX pcluster
read_cdfpart_cluster(int nc_file, const char *prefix);
#endif

#endif
