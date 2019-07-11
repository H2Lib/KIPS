#include "coulomb.h"
#include "spatialgeometry.h"
#include <stdio.h>
#include <stdlib.h>
#include "parameters.h"

#define M_NaCl -1.747564594633
#define a_NaCl 5.64

int 
main (int argc, char **argv) {
  preal bmin, bmax, x0, *x;
  uint dim, i, *j, k, *n, l, m, s, maxdepth, points, kernel;
  real eta, maxdiam, alpha, R, a, M, norm, error, t_setup;
  pspatialgeometry sg;
  pspatialcluster sroot;
  pblock broot;
  pclusterbasis cb;
  pparameters p;
  split sp;
  bool fw;
  pkernelmatrix km;
  pamatrix G;
  pavector q, y;
  size_t sz;
  pstopwatch sw;
  ph2matrix Gh2;
  
  init_kips (&argc, &argv);
  dim = 3;
  maxdepth = 8;
  maxdiam = 0.01;
  a = a_NaCl;
  alpha = 0.0;
  m = 3;
  points = 8000;
  eta = 1.0;
  
  kernel = askforint ("Splitting function? 0) SP, 1) SF, 2) SP2, 3) SP3, 4) DSP, 5) DSF", "", 3);
  if (kernel == 0) {
    sp = &split_SP;
  }
  else if (kernel == 1) {
    sp = &split_SF;
  }
  else if (kernel == 2) {
    sp = &split_SP2;
  }
  else if (kernel == 3) {
    sp = &split_SP3;
  }
  else if (kernel == 4) {
    sp = &split_DSP;
    alpha = askforreal ("Damping parameter?", "", 0.8/a);
  }
  else {
    sp = &split_DSF;
    alpha = askforreal ("Damping parameter?", "", 0.8/a);
  }
  R = askforreal ("Cutoff radius?", "", 9.0);
  
  sw = new_stopwatch();
  bmin = allocreal (dim);
  bmax = allocreal (dim);
  j = allocuint (dim);
  n = allocuint (dim);
  q = new_avector (points);
  y = new_avector (points);
  x = (preal *) allocmem (sizeof (preal) * points);
  x[0] = x0 = allocreal (points * dim);
  for (i=1; i<points; i++) {
    x0 += dim;
    x[i] = x0;
  }
  
  printf ("Creating bounding box and point charges\n");
  l = 0;
  s = 0;
  for (i=0; i<dim; i++) {
    n[i] = 20;
    bmin[i] = 0.0;
    bmax[i] = 10.0*a;
    q->v[0] = 1.0;
    // q->v[1] = -1.0;
    j[i] = 0;
  }
  l+=1;
  fw = true;
  while (fw) {
    for (i=0; i<dim; i++) {
      j[i]++;
      s++;
      if (j[i] < n[i]) {
        q->v[l] = (s%2 == 0) ? 1.0 : -1.0;
        l+= 1;
        fw = true;
        break;
      }
      else {
        j[i] = 0;
        s -= n[i];
        fw = false;
      }
    }
  }
  
  printf("Creating spatial cluster tree and kernelmatrix object\n");
  sg = new_spatialgeometry (dim, bmin, bmax);
  sroot = init_spatialgeometry (maxdepth, maxdiam, sg);
  printf ("%u spatial clusters\n", sroot->desc);
  p = setparameters_coulomb (sp, alpha, R, sg);
  km = new_kernelmatrix (dim, 0, m);
  km->g = &kernel_coulomb;
  km->data = p;
  
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

  printf("Creating raw H^2-matrix\n");
  start_stopwatch(sw);
  Gh2 = build_fromblock_h2matrix(broot, cb, cb);
  fill_h2matrix_kernelmatrix(km, Gh2);
  t_setup = stop_stopwatch(sw);
  sz = getsize_h2matrix(Gh2);
  printf("  %.2f seconds\n", t_setup);
  printf("  %.1f MB\n" "  %.1f KB/DoF\n", sz / 1048576.0, sz / 1024.0 / points);
  
  printf ("Creating locations for point charges\n");
  l = 0;
  s = 0;
  for (i=0; i<dim; i++) {
    n[i] = 20;
    x[0][i] = 0.0;
    j[i] = 0;
  }
  l+=1;
  fw = true;
  while (fw) {
    for (i=0; i<dim; i++) {
      j[i]++;
      s++;
      if (j[i] < n[i]) {
        for (k=0; k<dim; k++) {
          x[l][k] = j[k] * a * 0.5 ;
        }
        l+= 1;
        fw = true;
        break;
      }
      else {
        j[i] = 0;
        s -= n[i];
        fw = false;
      }
    }
  }
  
  printf ("Distributing point charges to clusters\n");
  initPoints_spatialgeometry (points, (pcreal *) x, sg);
  
  printf ("Updating H2-Matrix\n");
  start_stopwatch(sw);
  update_kernelmatrix (points, x, km);
  update_h2matrix_kernelmatrix (km, Gh2);
  t_setup = stop_stopwatch (sw);
  printf ("  %.6f seconds\n", t_setup);
  sz = getsize_h2matrix(Gh2);
  printf("  %.1f MB\n" "  %.1f KB/DoF\n", sz / 1048576.0, sz / 1024.0 / points);
  
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
  
  printf ("Computing Madelung constant of rock salt...\n");
  clear_avector (y);
  mvm_amatrix (1.0, false, G, q, y);
  M = a/points * (0.5*dotprod_avector (q,y) + selfPotential_coulomb(sp,alpha,q,R));
  printf ("\t via reference matrix: M = %10.8f    relative error: %10.8f\n", M, REAL_ABS((M-M_NaCl)/M_NaCl)); 
  clear_avector (y);
  mvm_h2matrix (1.0, false, Gh2, q, y);
  M = a/points * (0.5*dotprod_avector (q,y) + selfPotential_coulomb(sp,alpha,q,R));
  printf ("\t via compressed matrix: M = %10.8f    relative error: %10.8f\n", M, REAL_ABS((M-M_NaCl)/M_NaCl));
  
  printf ("Cleaning up\n");
  del_spatialgeometry (sg);
  del_spatialcluster (sroot);
  del_block (broot);
  del_clusterbasis (cb);
  freemem (p);
  del_kernelmatrix (km);
  freemem (j);
  freemem (n);
  del_avector (q);
  del_avector (y);
  del_amatrix (G);
  del_h2matrix (Gh2);
  
  return EXIT_SUCCESS;
}
