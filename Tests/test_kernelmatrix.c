
#include "spatialgeometry.h"
#include "kernelmatrix.h"
#include "parameters.h"
#include "basic.h"

#include <stdio.h>
#include <stdlib.h>

static real
kernel_newton(const real *xx, const real *yy, void *data)
{
  real norm2;

  (void) data;

  norm2 = REAL_SQR(xx[0] - yy[0]) + REAL_SQR(xx[1] - yy[1]) + REAL_SQR(xx[2] - yy[2]);

  return (norm2 == 0.0 ? 0.0 : 1.0 / REAL_SQRT(norm2));
}

static real
kernel_exp(const real *xx, const real *yy, void *data)
{
  real norm2;

  (void) data;

  norm2 = REAL_SQR(xx[0] - yy[0]) + REAL_SQR(xx[1] - yy[1]) + REAL_SQR(xx[2] - yy[2]);

  return REAL_EXP(-norm2);
}

static real
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
  real eta, norm, error, maxdiam, t_setup, *bmin, *bmax;
  uint i, j, dim, points, m, maxdepth;
  pspatialgeometry sg;
  pspatialcluster sroot;
  pblock broot;
  pclusterbasis cb;
  ph2matrix Gh2;
  pamatrix G;
  char kernel;
  pkernelmatrix km;
  pstopwatch sw;
  size_t sz;
  
  init_kips (&argc, &argv);
  
  sw = new_stopwatch();
  
  dim = askforint ("Spatial dimension?", "", 3);
  maxdepth = askforint ("Maximal depth of spatial cluster tree?", "", 6);
  maxdiam = askforreal ("Maximal diameter of leaf clusters?", "", 0.01);
  kernel = askforchar("Kernel function? N)ewton, L)ogarithmic, or E)xponential?", 
                      "h2lib_kernelfunc", "nle", 'n');
  points = askforint("Number of points?", "h2lib_kernelpoints", 2048);
  m = askforint("Interpolation order?", "h2lib_interorder", 3);

  eta = 1.0;
  
  printf("Creating kernelmatrix object for %u points, order %u\n", points, m);
  
  km = new_kernelmatrix(3, points, m);
  switch(kernel) {
  case 'e':
    printf("Exponential kernel function\n");
    km->g = &kernel_exp;
    break;
  case 'n':
    printf("Newton kernel function\n");
    km->g = &kernel_newton;
    break;
  default:
    printf("Logarithmic kernel function\n");
    km->g = &kernel_log;
  }
  
  printf ("Creating bounding box and random points\n");
  bmin = allocreal (dim);
  bmax = allocreal (dim);
  for (i=0; i<dim; i++) {
    bmin[i] = -1.0;
    bmax[i] = 1.0;
    for(j=0; j<points; j++) {
        km->x[j][i] = FIELD_RAND();
    }
  }
  
  printf("Creating spatial cluster tree\n");
  sg = new_spatialgeometry (dim, bmin, bmax);
  sroot = init_spatialgeometry (maxdepth, maxdiam, sg);
  printf ("%u spatial clusters\n", sroot->desc);
  
  printf ("Distributing random points to clusters\n");
  start_stopwatch (sw);
  initPoints_spatialgeometry (points, (const real **)km->x, sg);
  t_setup = stop_stopwatch (sw);
  printf ("%.6f seconds\n", t_setup);
  
  printf("Creating block tree\n");
  broot = buildh2periodic_block(sroot, sroot, eta, bmin, bmax);
  printf("%u blocks, depth %u\n",
		broot->desc, getdepth_block(broot));

  printf("Creating cluster basis\n");
  cb = build_fromcluster_clusterbasis(sroot);

  printf("Filling cluster basis\n");
  start_stopwatch(sw);
  fill_clusterbasis_kernelmatrix(km, cb);
  t_setup = stop_stopwatch(sw);
  sz = getsize_clusterbasis(cb);
  printf("  %.2f seconds\n" "  %.1f MB\n" "  %.1f KB/DoF\n", 
         t_setup, sz / 1048576.0, sz / 1024.0 / points);

  printf("Creating H^2-matrix\n");
  Gh2 = build_fromblock_h2matrix(broot, cb, cb);
  sz = getsize_h2matrix(Gh2);
  printf("  %.1f MB\n" "  %.1f KB/DoF\n", sz / 1048576.0, sz / 1024.0 / points);

  printf("Filling H^2-matrix\n");
  start_stopwatch(sw);
  fill_h2matrix_kernelmatrix(km, Gh2);
  t_setup = stop_stopwatch(sw);
  printf("  %.2f seconds\n", t_setup);

  printf("Computing norm\n");
  norm = norm2_h2matrix(Gh2);
  printf("  Spectral norm %.3e\n", norm);
  
  printf("Filling reference matrix\n");
  G = new_amatrix(points, points);
  start_stopwatch(sw);
  fillN_kernelmatrix(0, 0, km, G);
  t_setup = stop_stopwatch(sw);
  sz = getsize_amatrix(G);
  printf("  %.2f seconds\n" "  %.1f MB\n" "  %.1f KB/DoF\n", 
         t_setup, sz / 1048576.0, sz / 1024.0 / points);

  printf("Computing norm\n");
  norm = norm2_amatrix(G);
  printf("  Spectral norm %.3e\n", norm);
  
  printf("Computing approximation error\n");
  error = norm2diff_amatrix_h2matrix(Gh2, G);
  printf("  Spectral error %.3e (%.3e)\n", error, error/norm);
  
  printf ("Cleaning up\n");
  del_amatrix (G);
  del_h2matrix (Gh2);
  del_clusterbasis (cb);
  del_kernelmatrix (km);
  del_block (broot);
  del_spatialgeometry (sg);
  del_spatialcluster (sroot);
  printf ("%u spatial clusters still active\n", 
          getactives_spatialcluster());
  printf ("%u blocks still active\n", getactives_block());
  
  uninit_kips();

  return EXIT_SUCCESS;
}
