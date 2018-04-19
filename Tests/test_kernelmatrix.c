
#include "kernelmatrix.h"

#include <stdio.h>

#include "basic.h"
#include "h2matrix.h"
#include "parameters.h"
#include "amatrix.h"
#include "eigensolvers.h"

static field
kernel_newton(const real *xx, const real *yy, void *data)
{
  real norm2;

  (void) data;

  norm2 = REAL_SQR(xx[0] - yy[0]) + REAL_SQR(xx[1] - yy[1]) + REAL_SQR(xx[2] - yy[2]);

  return (norm2 == 0.0 ? 0.0 : 1.0 / REAL_SQRT(norm2));
}

static field
kernel_exp(const real *xx, const real *yy, void *data)
{
  real norm2;

  (void) data;

  norm2 = REAL_SQR(xx[0] - yy[0]) + REAL_SQR(xx[1] - yy[1]) + REAL_SQR(xx[2] - yy[2]);

  return REAL_EXP(-norm2);
}

static field
kernel_log(const real *xx, const real *yy, void *data)
{
  real norm2;

  (void) data;

  norm2 = REAL_SQR(xx[0] - yy[0]) + REAL_SQR(xx[1] - yy[1]) + REAL_SQR(xx[2] - yy[2]);

  return (norm2 == 0.0 ? 0.0 : -0.5*REAL_LOG(norm2));
}

int
main(int argc, char **argv)
{
  pkernelmatrix km;
  pclustergeometry cg;
  pcluster root;
  pblock broot;
  pclusterbasis cb;
  ph2matrix Gh1;
  pamatrix G;
  pstopwatch sw;
  char kernel;
  uint points;
  uint m, leafsize;
  uint *idx;
  size_t sz;
  real eta, mindiam;
  real t_setup, norm, error;
  uint i;
  
  init_kips(&argc, &argv);

  sw = new_stopwatch();

  kernel = askforchar("Kernel function? N)ewton, L)ogarithmic, or E)xponential?", "h2lib_kernelfunc", "nle", 'n');

  points = askforint("Number of points?", "h2lib_kernelpoints", 2048);

  m = askforint("Interpolation order?", "h2lib_interorder", 3);

  leafsize = askforint("Cluster resolution?", "h2lib_leafsize", 2*m*m);

  eta = 2.0;
  
  (void) printf("Creating kernelmatrix object for %u points, order %u\n",
		points, m);
  km = new_kernelmatrix(3, points, m);
  switch(kernel) {
  case 'e':
    (void) printf("  Exponential kernel function\n");
    km->kernel = kernel_exp;
    break;
  case 'n':
    (void) printf("  Newton kernel function\n");
    km->kernel = kernel_newton;
    break;
  default:
    (void) printf("  Logarithmic kernel function\n");
    km->kernel = kernel_log;
  }
  for(i=0; i<points; i++) {
    km->x[i][0] = FIELD_RAND();	/* Random points in [-1,1]^2 */
    km->x[i][1] = FIELD_RAND();
    km->x[i][2] = FIELD_RAND();
  }

  (void) printf("Creating clustergeometry object\n");
  cg = creategeometry_kernelmatrix(km);

  (void) printf("Creating cluster tree\n");
  idx = (uint *) allocmem(sizeof(uint) * points);
  for(i=0; i<points; i++)
    idx[i] = i;
  root = buildgeometric_clustergeometry(cg, 100, leafsize, points, idx);
  (void) printf("  %u clusters, depth %u\n",
		root->desc, getdepth_cluster(root));

  (void) printf("Checking for degenerate clusters\n");
  mindiam = getmindiam_cluster(root);
  (void) printf("  Minimal diameter %.4e\n", mindiam);
  if(mindiam == 0.0) {
    (void) printf("  That's too small, please try a larger cluster resolution\n");
    return 1;
  }
  
  (void) printf("Creating block tree\n");
  broot = buildh2std_block(root, root, eta);
  (void) printf("  %u blocks, depth %u\n",
		broot->desc, getdepth_block(broot));

  (void) printf("Creating cluster basis\n");
  cb = build_fromcluster_clusterbasis(root);

  (void) printf("Filling cluster basis\n");
  start_stopwatch(sw);
  fill_clusterbasis_kernelmatrix(km, cb);
  t_setup = stop_stopwatch(sw);
  sz = getsize_clusterbasis(cb);
  (void) printf("  %.2f seconds\n"
		"  %.1f MB\n"
		"  %.1f KB/DoF\n",
		t_setup, sz / 1048576.0, sz / 1024.0 / points);

  (void) printf("Creating H^2-matrix\n");
  Gh1 = build_fromblock_h2matrix(broot, cb, cb);
  sz = getsize_h2matrix(Gh1);
  (void) printf("  %.1f MB\n"
		"  %.1f KB/DoF\n",
		sz / 1048576.0, sz / 1024.0 / points);

  (void) printf("Filling H^2-matrix\n");
  start_stopwatch(sw);
  fill_h2matrix_kernelmatrix(km, Gh1);
  t_setup = stop_stopwatch(sw);
  (void) printf("  %.2f seconds\n",
		t_setup);

  (void) printf("Computing norm\n");
  norm = norm2_h2matrix(Gh1);
  (void) printf("  Spectral norm %.3e\n",
		norm);
  
  (void) printf("Filling reference matrix\n");
  G = new_amatrix(points, points);
  start_stopwatch(sw);
  fillN_kernelmatrix(0, 0, km, G);
  t_setup = stop_stopwatch(sw);
  sz = getsize_amatrix(G);
  (void) printf("  %.2f seconds\n"
		"  %.1f MB\n"
		"  %.1f KB/DoF\n",
		t_setup, sz / 1048576.0, sz / 1024.0 / points);

  (void) printf("Computing norm\n");
  norm = norm2_amatrix(G);
  (void) printf("  Spectral norm %.3e\n",
		norm);
  
  (void) printf("Computing approximation error\n");
  error = norm2diff_amatrix_h2matrix(Gh1, G);
  (void) printf("  Spectral error %.3e (%.3e)\n",
		error, error/norm);
  
  del_amatrix (G);
  del_h2matrix (Gh1);
  del_clustergeometry (cg);
  del_clusterbasis (cb);
  del_kernelmatrix (km);
  del_block (broot);
  del_cluster (root);
  free (idx);
  
  (void) printf("----------------------------------------\n"
		"  %u matrices and\n"
		"  %u vectors and\n"
        "  %u clusters and\n"
        "  %u blocks still active\n", 
		getactives_amatrix(), getactives_avector(), 
        getactives_cluster(), getactives_block());

  uninit_kips();

  return 0;
}
