
/* ------------------------------------------------------------
   This is the file "cluster.c" of the H2Lib package.
   All rights reserved, Steffen Boerm 2010
   ------------------------------------------------------------ */

#include "cluster.h"

#include "basic.h"

#include <assert.h>
#include <math.h>

#ifdef USE_NETCDF
#include <stdio.h>
#include <string.h>
#include <netcdf.h>
#endif

static uint active_cluster = 0;

/* ------------------------------------------------------------
 * Constructors and destructors
 * ------------------------------------------------------------ */

pcluster
new_cluster(uint size, uint *idx, uint sons, uint dim)
{
  pcluster t;
  uint i;

  t = (pcluster) allocmem(sizeof(cluster));
  t->size = size;
  t->idx = idx;
  t->sons = sons;
  t->dim = dim;
  t->desc = 0;

  t->son = NULL;
  if(sons > 0) {
    t->son = (pcluster *) allocmem((size_t) sizeof(cluster) * sons);
    for(i=0; i<sons; i++)
      t->son[i] = NULL;
  }

  t->bmin = allocreal(dim);
  t->bmax = allocreal(dim);

#ifdef USE_OPENMP
#pragma omp atomic
#endif
  active_cluster++;

  return t;
}

void
update_cluster(pcluster t)
{
  uint desc;
  uint i;

  desc = 1;
  for(i=0; i<t->sons; i++)
    desc += t->son[i]->desc;

  t->desc = desc;
}

void
del_cluster(pcluster t)
{
  uint i;
      
  freemem(t->bmin);
  freemem(t->bmax);

  if(t->sons > 0) {
    for(i=0; i<t->sons; i++)
      del_cluster(t->son[i]);
    freemem(t->son);
  }

  freemem(t);

  assert(active_cluster > 0);

#ifdef USE_OPENMP
#pragma omp atomic
#endif
  active_cluster--;
}

void
setsons_cluster(pcluster t, uint sons)
{
  uint i;

  assert(t->sons == 0);

  t->sons = sons;
  t->son = (pcluster *) allocmem(sizeof(pcluster) * sons);

  for(i=0; i<sons; i++)
    t->son[i] = 0;
}

/* ------------------------------------------------------------
 * Bounding box
 * ------------------------------------------------------------ */

real
diam_cluster(pccluster t)
{
  real diam2;
  uint i;

  diam2 = 0.0;
  for(i=0; i<t->dim; i++)
    diam2 += REAL_SQR(t->bmax[i] - t->bmin[i]);

  return REAL_SQRT(diam2);
}

real
dist_cluster(pccluster t, pccluster s)
{
  real dist2;
  uint i;

  assert(t->dim == s->dim);

  dist2 = 0.0;
  for(i=0; i<t->dim; i++)
    if(t->bmax[i] < s->bmin[i])
      dist2 += REAL_SQR(s->bmin[i] - t->bmax[i]);
    else if(s->bmax[i] < t->bmin[i])
      dist2 += REAL_SQR(t->bmin[i] - s->bmax[i]);

  return REAL_SQRT(dist2);
}

/* ------------------------------------------------------------
 * Statistics
 * ------------------------------------------------------------ */

uint
getactives_cluster()
{
  return active_cluster;
}

size_t
getsize_cluster(pccluster t)
{
  size_t sz;
  uint i;

  sz = (size_t) sizeof(cluster) + 2 * t->dim * sizeof(real);
  if(t->sons > 0) {
    sz += (size_t) sizeof(pcluster) * t->sons;
    for(i=0; i<t->sons; i++)
      sz += getsize_cluster(t->son[i]);
  }

  return sz;
}

uint
getdepth_cluster(pccluster t)
{
  uint d, d1;
  uint i;

  d = 0;

  if(t->sons) {
    for(i=0; i<t->sons; i++) {
      d1 = getdepth_cluster(t->son[i]);
      if(d1 > d)
	d = d1;
    }
    d++;
  }
  
  return d;
}

/* ------------------------------------------------------------
 * Hierarchical iterator
 * ------------------------------------------------------------ */

void
iterate_cluster(pccluster t, uint tname,
		void (*pre)(pccluster t, uint tname, void *data),
		void (*post)(pccluster t, uint tname, void *data),
		void *data)
{
#ifdef USE_OPENMP
  iterate_parallel_cluster(t, tname, max_pardepth, pre, post, data);
#else
  uint tname1;
  uint i;

  if(pre)
    pre(t, tname, data);

  tname1 = tname + 1;
  for(i=0; i<t->sons; i++) {
    iterate_cluster(t->son[i], tname1,
		    pre, post, data);

    tname1 += t->son[i]->desc;
  }
  assert(tname1 == tname + t->desc);
  
  if(post)
    post(t, tname, data);
#endif
}

void
iterate_parallel_cluster(pccluster t, uint tname, uint pardepth,
			 void (*pre)(pccluster t, uint tname, void *data),
			 void (*post)(pccluster t, uint tname, void *data),
			 void *data)
{
  uint *tname1, tname2;
#ifdef USE_OPENMP
  uint nthreads;		/* HACK: Solaris workaround */
#endif
  uint i;
  int p;

  if(pre)
    pre(t, tname, data);

  if(t->sons > 0) {
    tname1 = (uint *) allocmem((size_t) sizeof(uint) * t->sons);

    tname2 = tname + 1;
    for(i=0; i<t->sons; i++) {
      tname1[i] = tname2;
      tname2 += t->son[i]->desc;
    }
    assert(tname2 == tname + t->desc);

#ifdef USE_OPENMP
    nthreads = t->sons;
    (void) nthreads;
#pragma omp parallel for if(pardepth > 0), num_threads(nthreads)
#endif
    for(p=0; p<(int) t->sons; p++)
      iterate_parallel_cluster(t->son[p], tname1[p],
			       (pardepth > 0 ? pardepth-1 : 0),
			       pre, post, data);

    freemem(tname1);
  }
  
  if(post)
    post(t, tname, data);
}

/* ------------------------------------------------------------
 * Enumeration
 * ------------------------------------------------------------ */

static void
enumerate(pcluster t, uint tname, pcluster *tn)
{
  uint tname1;
  uint i;

  tn[tname] = t;

  tname1 = tname + 1;
  for(i=0; i<t->sons; i++) {
    enumerate(t->son[i], tname1, tn);

    tname1 += t->son[i]->desc;
  }
  assert(tname1 = tname + t->desc);
}

pcluster *
enumerate_cluster(pcluster t)
{
  pcluster *tn;
  
  tn = (pcluster *) allocmem((size_t) sizeof(pcluster) * t->desc);

  enumerate(t, 0, tn);

  return tn;
}

/* ------------------------------------------------------------
 * File I/O
 * ------------------------------------------------------------ */

#ifdef USE_NETCDF
static void
write_count(pccluster t, size_t *clusters, size_t *coeffs)
{
  uint i;

  /* Increase cluster counter */
  (*clusters)++;

  /* Add number of coefficients of transfer matrix */
  (*coeffs) += 2 * t->dim;

  /* Handle sons */
  for(i=0; i<t->sons; i++)
    write_count(t->son[i], clusters, coeffs);
}

static void
write_cdf(pccluster t,
	  size_t clusters, size_t coeffs,
	  size_t *clusteridx, size_t *coeffidx,
	  int nc_file,
	  int nc_sons, int nc_size, int nc_coeff)
{
  size_t start, count;
  ptrdiff_t stride;
  int val, result;
  uint i;

  assert(*clusteridx <= clusters);

  /* Write number of sons to nc_sons[*clusteridx] */
  start = *clusteridx;
  count = 1;
  stride = 1;
  val = t->sons;
  result = nc_put_vars(nc_file, nc_sons, &start, &count, &stride, &val);
  assert(result == NC_NOERR);

  /* Write size of cluster to nc_size[*clusteridx] */
  val = t->size;
  result = nc_put_vars(nc_file, nc_size, &start, &count, &stride, &val);
  assert(result == NC_NOERR);

  /* Increase cluster index */
  (*clusteridx)++;

  /* Handle sons */
  for(i=0; i<t->sons; i++)
    write_cdf(t->son[i], clusters, coeffs, clusteridx, coeffidx,
	      nc_file, nc_sons, nc_size, nc_coeff);

  /* Write bounding box */
  start = *coeffidx;
  assert(start + 2*t->dim <= coeffs);
  count = t->dim;
  result = nc_put_vars(nc_file, nc_coeff, &start, &count, &stride, t->bmin);
  assert(result == NC_NOERR);
  start += t->dim;

  result = nc_put_vars(nc_file, nc_coeff, &start, &count, &stride, t->bmax);
  assert(result == NC_NOERR);
  start += t->dim;
  (*coeffidx) = start;
}

void
write_cdf_cluster(pccluster t, const char *name)
{
  size_t clusters, clusteridx;
  size_t coeffs, coeffidx;
  int nc_file, nc_sons, nc_size, nc_idx, nc_coeff;
  int nc_clusters, nc_coeffs, nc_totalsize, nc_dim;
  int result;

  /* Count number of clusters and coefficients */
  clusters = 0;
  coeffs = 0;
  write_count(t, &clusters, &coeffs);

  /* Create NetCDF file */
  result = nc_create(name, NC_64BIT_OFFSET, &nc_file);
  assert(result == NC_NOERR);

  /* Define "clusters" dimension */
  result = nc_def_dim(nc_file, "clusters", clusters, &nc_clusters);
  assert(result == NC_NOERR);

  /* Define "coeffs" dimension */
  result = nc_def_dim(nc_file, "coeffs", coeffs, &nc_coeffs);
  assert(result == NC_NOERR);

  /* Define "totalsize" dimension */
  result = nc_def_dim(nc_file, "totalsize", t->size, &nc_totalsize);
  assert(result == NC_NOERR);

  /* Define "dim" dimension */
  result = nc_def_dim(nc_file, "dim", t->dim, &nc_dim);
  assert(result == NC_NOERR);

  /* Define "sons" variable */
  result = nc_def_var(nc_file, "sons", NC_INT, 1, &nc_clusters, &nc_sons);
  assert(result == NC_NOERR);

  /* Define "size" variable */
  result = nc_def_var(nc_file, "size", NC_INT, 1, &nc_clusters, &nc_size);
  assert(result == NC_NOERR);

  /* Define "idx" variable */
  result = nc_def_var(nc_file, "idx", NC_INT, 1, &nc_totalsize, &nc_idx);
  assert(result == NC_NOERR);

  /* Define "coeff" variable */
  result = nc_def_var(nc_file, "coeff", NC_DOUBLE, 1, &nc_coeffs, &nc_coeff);
  assert(result == NC_NOERR);

  /* Finish NetCDF define mode */
  result = nc_enddef(nc_file);
  assert(result == NC_NOERR);

  /* Write index to NetCDF variable */
  result = nc_put_var(nc_file, nc_idx, t->idx);

  /* Write coefficiens to NetCDF variables */
  clusteridx = 0;
  coeffidx = 0;
  write_cdf(t, clusters, coeffs, &clusteridx, &coeffidx,
	    nc_file, nc_sons, nc_size, nc_coeff);
  assert(clusteridx == clusters);
  assert(coeffidx == coeffs);

  /* Close file */
  result = nc_close(nc_file);
  assert(result == NC_NOERR);
}

static void
prefix_name(char *buf, int bufsize, const char *prefix, const char *name)
{
  if(prefix)
    snprintf(buf, bufsize, "%s_%s", prefix, name);
  else
    strncpy(buf, name, bufsize);
}

void
write_cdfpart_cluster(pccluster t, int nc_file, const char *prefix)
{
  size_t clusters, clusteridx;
  size_t coeffs, coeffidx;
  char *buf;
  int bufsize;
  int nc_sons, nc_size, nc_idx, nc_coeff;
  int nc_clusters, nc_coeffs, nc_totalsize, nc_dim;
  int result;

  /* Prepare buffer for prefixed names */
  bufsize = strlen(prefix) + 16;
  buf = (char *) allocmem(sizeof(char) * bufsize);

  /* Count number of clusters and coefficients */
  clusters = 0;
  coeffs = 0;
  write_count(t, &clusters, &coeffs);

  /* Switch NetCDF file to define mode */
  result = nc_redef(nc_file);
  assert(result == NC_NOERR || result == NC_EINDEFINE);

  /* Define "clusters" dimension */
  prefix_name(buf, bufsize, prefix, "clusters");
  result = nc_def_dim(nc_file, buf, clusters, &nc_clusters);
  assert(result == NC_NOERR);

  /* Define "coeffs" dimension */
  prefix_name(buf, bufsize, prefix, "coeffs");
  result = nc_def_dim(nc_file, buf, coeffs, &nc_coeffs);
  assert(result == NC_NOERR);

  /* Define "totalsize" dimension */
  prefix_name(buf, bufsize, prefix, "totalsize");
  result = nc_def_dim(nc_file, buf, t->size, &nc_totalsize);
  assert(result == NC_NOERR);

  /* Define "dim" dimension */
  prefix_name(buf, bufsize, prefix, "dim");
  result = nc_def_dim(nc_file, buf, t->dim, &nc_dim);
  assert(result == NC_NOERR);

  /* Define "sons" variable */
  prefix_name(buf, bufsize, prefix, "sons");
  result = nc_def_var(nc_file, buf, NC_INT, 1, &nc_clusters, &nc_sons);
  assert(result == NC_NOERR);

  /* Define "size" variable */
  prefix_name(buf, bufsize, prefix, "size");
  result = nc_def_var(nc_file, buf, NC_INT, 1, &nc_clusters, &nc_size);
  assert(result == NC_NOERR);

  /* Define "idx" variable */
  prefix_name(buf, bufsize, prefix, "idx");
  result = nc_def_var(nc_file, buf, NC_INT, 1, &nc_totalsize, &nc_idx);
  assert(result == NC_NOERR);

  /* Define "coeff" variable */
  prefix_name(buf, bufsize, prefix, "coeff");
  result = nc_def_var(nc_file, buf, NC_DOUBLE, 1, &nc_coeffs, &nc_coeff);
  assert(result == NC_NOERR);

  /* Finish NetCDF define mode */
  result = nc_enddef(nc_file);
  assert(result == NC_NOERR);

  /* Write index to NetCDF variable */
  result = nc_put_var(nc_file, nc_idx, t->idx);

  /* Write coefficiencs to NetCDF variables */
  clusteridx = 0;
  coeffidx = 0;
  write_cdf(t, clusters, coeffs, &clusteridx, &coeffidx,
	    nc_file, nc_sons, nc_size, nc_coeff);
  assert(clusteridx == clusters);
  assert(coeffidx == coeffs);

  /* Clean up */
  nc_sync(nc_file);
  freemem(buf);
}

static pcluster
read_cdf_part(int nc_file, size_t clusters, size_t coeffs,
	      int nc_sons, int nc_size, int nc_coeff,
	      uint *idx, int dim,
	      size_t *clusteridx, size_t *coeffidx)
{
  pcluster t, t1;
  uint *idx1;
  uint size;
  uint sons;
  uint i;
  size_t start, count;
  ptrdiff_t stride;
  int val, result;

  /* Get number of sons */
  start = *clusteridx;
  count = 1;
  stride = 1;
  result = nc_get_vars(nc_file, nc_sons, &start, &count, &stride, &val);
  assert(result == NC_NOERR);
  sons = val;

  /* Get size of cluster */
  result = nc_get_vars(nc_file, nc_size, &start, &count, &stride, &val);
  assert(result == NC_NOERR);
  size = val;

  /* Create new cluster */
  t = new_cluster(size, idx, sons, dim);

  /* Increase cluster index */
  (*clusteridx)++;

  /* Handle sons */
  if(sons > 0) {
    idx1 = idx;
    for(i=0; i<sons; i++) {
      t1 = read_cdf_part(nc_file, clusters, coeffs,
			 nc_sons, nc_size, nc_coeff,
			 idx1, dim, clusteridx, coeffidx);
      t->son[i] = t1;
      
      idx1 += t1->size;
    }
    assert(idx1 == idx + size);
  }

  /* Get bounding box */
  start = (*coeffidx);
  count = dim;
  result = nc_get_vars(nc_file, nc_coeff, &start, &count, &stride, t->bmin);
  start += dim;

  result = nc_get_vars(nc_file, nc_coeff, &start, &count, &stride, t->bmax);
  start += dim;
  (*coeffidx) = start;

  /* Finish initialization */
  update_cluster(t);

  return t;
}

pcluster
read_cdf_cluster(const char *name)
{
  pcluster t;
  uint *idx;
  size_t clusters, clusteridx;
  size_t coeffs, coeffidx;
  char dimname[NC_MAX_NAME+1];
  int nc_file, nc_sons, nc_size, nc_idx, nc_coeff;
  int nc_clusters, nc_coeffs, nc_totalsize, nc_dim;
  size_t dim, totalsize;
  int result;

  /* Open NetCDF file */
  result = nc_open(name, NC_NOWRITE, &nc_file);
  assert(result == NC_NOERR);

  /* Get "clusters" dimension */
  result = nc_inq_dimid(nc_file, "clusters", &nc_clusters);
  assert(result == NC_NOERR);
  result = nc_inq_dim(nc_file, nc_clusters, dimname, &clusters);
  assert(result == NC_NOERR);

  /* Get "coeffs" dimension */
  result = nc_inq_dimid(nc_file, "coeffs", &nc_coeffs);
  assert(result == NC_NOERR);
  result = nc_inq_dim(nc_file, nc_coeffs, dimname, &coeffs);
  assert(result == NC_NOERR);

  /* Get "totalsize" dimension */
  result = nc_inq_dimid(nc_file, "totalsize", &nc_totalsize);
  assert(result == NC_NOERR);
  result = nc_inq_dim(nc_file, nc_totalsize, dimname, &totalsize);
  assert(result == NC_NOERR);

  /* Get "dim" dimension */
  result = nc_inq_dimid(nc_file, "dim", &nc_dim);
  assert(result == NC_NOERR);
  result = nc_inq_dim(nc_file, nc_dim, dimname, &dim);
  assert(result == NC_NOERR);

  /* Get "sons" variable */
  result = nc_inq_varid(nc_file, "sons", &nc_sons);
  assert(result == NC_NOERR);

  /* Get "size" variable */
  result = nc_inq_varid(nc_file, "size", &nc_size);
  assert(result == NC_NOERR);

  /* Get "idx" variable */
  result = nc_inq_varid(nc_file, "idx", &nc_idx);
  assert(result == NC_NOERR);

  /* Get "coeff" variable */
  result = nc_inq_varid(nc_file, "coeff", &nc_coeff);
  assert(result == NC_NOERR);

  /* Read index */
  idx = (uint *) allocmem(sizeof(uint) * totalsize);
  nc_get_var(nc_file, nc_idx, idx);

  /* Read coefficients from NetCDF variables */
  clusteridx = 0;
  coeffidx = 0;
  t = read_cdf_part(nc_file, clusters, coeffs,
		    nc_sons, nc_size, nc_coeff,
		    idx, dim, &clusteridx, &coeffidx);
  assert(clusteridx == clusters);
  assert(coeffidx == coeffs);

  /* Close NetCDF file */
  nc_close(nc_file);

  return t;
}

pcluster
read_cdfpart_cluster(int nc_file, const char *prefix)
{
  pcluster t;
  uint *idx;
  size_t clusters, clusteridx;
  size_t coeffs, coeffidx;
  char dimname[NC_MAX_NAME+1];
  char *buf;
  int bufsize;
  int nc_sons, nc_size, nc_idx, nc_coeff;
  int nc_clusters, nc_coeffs, nc_totalsize, nc_dim;
  size_t dim, totalsize;
  int result;

  /* Prepare buffer for prefixed names */
  bufsize = strlen(prefix) + 16;
  buf = (char *) allocmem(sizeof(char) * bufsize);

  /* Get "clusters" dimension */
  prefix_name(buf, bufsize, prefix, "clusters");
  result = nc_inq_dimid(nc_file, buf, &nc_clusters);
  assert(result == NC_NOERR);
  result = nc_inq_dim(nc_file, nc_clusters, dimname, &clusters);
  assert(result == NC_NOERR);

  /* Get "coeffs" dimension */
  prefix_name(buf, bufsize, prefix, "coeffs");
  result = nc_inq_dimid(nc_file, buf, &nc_coeffs);
  assert(result == NC_NOERR);
  result = nc_inq_dim(nc_file, nc_coeffs, dimname, &coeffs);
  assert(result == NC_NOERR);

  /* Get "totalsize" dimension */
  prefix_name(buf, bufsize, prefix, "totalsize");
  result = nc_inq_dimid(nc_file, buf, &nc_totalsize);
  assert(result == NC_NOERR);
  result = nc_inq_dim(nc_file, nc_totalsize, dimname, &totalsize);
  assert(result == NC_NOERR);

  /* Get "dim" dimension */
  prefix_name(buf, bufsize, prefix, "dim");
  result = nc_inq_dimid(nc_file, buf, &nc_dim);
  assert(result == NC_NOERR);
  result = nc_inq_dim(nc_file, nc_dim, dimname, &dim);
  assert(result == NC_NOERR);

  /* Get "sons" variable */
  prefix_name(buf, bufsize, prefix, "sons");
  result = nc_inq_varid(nc_file, buf, &nc_sons);
  assert(result == NC_NOERR);

  /* Get "size" variable */
  prefix_name(buf, bufsize, prefix, "size");
  result = nc_inq_varid(nc_file, buf, &nc_size);
  assert(result == NC_NOERR);

  /* Get "idx" variable */
  prefix_name(buf, bufsize, prefix, "idx");
  result = nc_inq_varid(nc_file, buf, &nc_idx);
  assert(result == NC_NOERR);

  /* Get "coeff" variable */
  prefix_name(buf, bufsize, prefix, "coeff");
  result = nc_inq_varid(nc_file, buf, &nc_coeff);
  assert(result == NC_NOERR);

  /* Read index */
  idx = (uint *) allocmem(sizeof(uint) * totalsize);
  nc_get_var(nc_file, nc_idx, idx);

  /* Read coefficients from NetCDF variables */
  clusteridx = 0;
  coeffidx = 0;
  t = read_cdf_part(nc_file, clusters, coeffs,
		    nc_sons, nc_size, nc_coeff,
		    idx, dim, &clusteridx, &coeffidx);
  assert(clusteridx == clusters);
  assert(coeffidx == coeffs);

  /* Clean up */
  freemem(buf);

  return t;
}
#endif
