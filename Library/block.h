
/* ------------------------------------------------------------
 * This is the file "block.h" of the KIPS package.
 * All rights reserved, Steffen Boerm 2010
 * ------------------------------------------------------------ */

#ifndef BLOCK_H
#define BLOCK_H

typedef struct _admflag admflag;
typedef admflag *padmflag;

typedef struct _block block;
typedef block *pblock;
typedef const block *pcblock;

#include "spatialcluster.h"
#include "settings.h"

#include <stdlib.h>

#ifdef USE_CAIRO
#include <cairo.h>
#endif

struct _block
{
  pspatialcluster rc;			/* Row cluster */
  pspatialcluster cc;			/* Column cluster */

  bool adm;			/* Admissibility */

  uint rsons;			/* Number of row sons */
  uint csons;			/* Number of column sons */
  pblock *son;			/* Pointers to sons */

  uint desc;			/* Number of descendants */
};

/* ------------------------------------------------------------
 * Constructors and destructors
 * ------------------------------------------------------------ */

HEADER_PREFIX pblock
new_block(pspatialcluster rc, pspatialcluster cc, uint rsons, uint csons);

HEADER_PREFIX void
update_block(pblock b);

HEADER_PREFIX void
del_block(pblock b);

/* ------------------------------------------------------------
 * Statistics
 * ------------------------------------------------------------ */

HEADER_PREFIX uint
getactives_block();

HEADER_PREFIX size_t
getsize_block(pcblock t);

HEADER_PREFIX uint
getdepth_block(pcblock b);

/* ------------------------------------------------------------
 * Drawing
 * ------------------------------------------------------------ */

#ifdef USE_CAIRO
HEADER_PREFIX void
cairodraw_block(cairo_t *cr, pcblock b, int levels);
#endif

/* ------------------------------------------------------------
 * Standard admissibility conditions,
 * "data" points to the admissibility parameter eta
 * ------------------------------------------------------------ */

HEADER_PREFIX bool
h2std_admissibility(pspatialcluster rc, pspatialcluster cc, void *data);

HEADER_PREFIX bool
hstd_admissibility(pspatialcluster rc, pspatialcluster cc, void *data);

/* Standard admissibility condition for a periodic system. 
 * The parameter "data" is supposed to contain the admissibility 
 * parameter eta, the vector bmin of the minimal coordinates of the unit cell and the 
 * vector bmax of the corresponding maximal coordinates, in this order.*/
HEADER_PREFIX bool
h2periodic_admissibility (pspatialcluster rc, pspatialcluster cc, void *data);

/* ------------------------------------------------------------
 * Building standard block trees
 * ------------------------------------------------------------ */

HEADER_PREFIX pblock
build_block(pspatialcluster rc, pspatialcluster cc,
	    bool (*admissible)(pspatialcluster rc, pspatialcluster cc, void *data),
	    void *data,
	    uint levels, bool strict);

HEADER_PREFIX pblock
buildh2std_block(pspatialcluster rc, pspatialcluster cc, real eta);

HEADER_PREFIX pblock
buildhstd_block(pspatialcluster rc, pspatialcluster cc, real eta);

HEADER_PREFIX pblock
buildh2periodic_block (pspatialcluster rc, pspatialcluster cc, real eta, preal bmin, preal bmax);

#endif
