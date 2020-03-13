#include "coulomb.h"
#include "lj.h"
#include "tip4p.h"
#include "spatialgeometry.h"
#include <stdio.h>
#include <stdlib.h>
#include "parameters.h"
#include "quaternion.h"

int 
main (int argc, char **argv) {
  preal bmin, bmax, *xn, *xc, w, I, v;//, r, com;
  uint dim, i, j, m, maxdepth, n;
  real eta, mindiam, R, E, E_rot, E_trans, sig, eps, alpha, delta, time;
  pspatialgeometry sg;
  pspatialcluster sroot;
  pblock broot;
  pclusterbasis rbc, cbc, rblj, cblj;
  pcoulomb c;
  plj lj;
  split s, sd;
  pkernelmatrix kc, klj;
  ph2matrix Vc, Vlj, Fc, Flj;
  ptip4p t;
  pavector fc, flj;
  
  init_kips (&argc, &argv);
  
  dim = 3;
  maxdepth = 0;
  mindiam = RHH_TIP4P;
  alpha = 0.0;
  m = 3;
  eta = 1.0;
  sig = SIG_TIP4P;
  eps = EPS_TIP4P;
  s = &split_cutoff;
  sd = &derivative_cutoff;
  R = 5.0e-8;
  bmin = allocreal (dim);
  bmax = allocreal (dim);
  
  
  printf ("Setting up bounding box and cluster and block trees\n");
  for (i=0; i<dim; i++) {
    bmin[i] = -1.0e-7;
    bmax[i] = 1.0e-7;
  }
  
  sg = new_spatialgeometry (dim, bmin, bmax);
  sroot = init_spatialgeometry (maxdepth, mindiam, sg);
  c = setParametersGradient_coulomb (s, sd, alpha, R, sg);
  lj = setParameters_lj (R, sig, eps, sg);
  kc = new_kernelmatrix (dim, m);
  klj = new_kernelmatrix (dim, m);
  kc->g = &kernel_coulomb;
  kc->f = &gradient_coulomb;
  kc->data = c;
  klj->g = &kernel_lj;
  klj->f = &gradient_lj;
  klj->data = lj;
  
  broot = buildh2periodic_block(sroot, sroot, eta, bmin, bmax);

  cbc = build_fromcluster_clusterbasis (sroot);
  cblj = build_fromcluster_clusterbasis (sroot);
  rbc = build_fromcluster_clusterbasis (sroot);
  rblj = build_fromcluster_clusterbasis (sroot);

  printf ("Setting up H2-matrix structures\n");

  fill_clusterbasis_kernelmatrix (false, kc, cbc);
  fill_clusterbasis_kernelmatrix (false, klj, cblj);
  fill_clusterbasis_kernelmatrix (true, kc, rbc);
  fill_clusterbasis_kernelmatrix (true, klj, rblj);

  Vc = build_fromblock_h2matrix (broot, cbc, cbc);
  Vlj = build_fromblock_h2matrix (broot, cblj, cblj);
  Fc = build_fromblock_h2matrix (broot, rbc, cbc);
  Flj = build_fromblock_h2matrix (broot, rblj, cblj);
  fill_h2matrix_kernelmatrix (false, kc, Vc);
  fill_h2matrix_kernelmatrix (false, klj, Vlj);
  fill_h2matrix_kernelmatrix (true, kc, Fc);
  fill_h2matrix_kernelmatrix (true, klj, Flj);
  
  printf ("Reading input files\n");
  
  t = inputfile_tip4p ("./Experiments/coordinates.tip4p");
  n = t->n;
  fc = new_avector (3 * n * 3);
  flj = new_avector (n * 3);
  update_kernelmatrix (3 * n, t->xc, t->xcmol, kc);
  update_kernelmatrix (n, t->xn, t->xnmol, klj);
  fullForce_coulomb (t->q, kc, sg, Fc, fc);
  fullForce_lj (klj, sg, Flj, flj);
  calcForce_tip4p (fc, flj, t);
  calcTorque_tip4p (fc, flj, t);
  initMolecules_tip4p (t);
  
  time = 0.0;
  E = initVelocitiesFromFile_molecule ("./Experiments/velocities.tip4p", n, t->mol);
  printf ("t = %lf  T = %lf\n", time, 2.0 * E / ((6 * n - 3) * K_B));
  
  delta = 1.0e-14;
  for (i=0; i<10; i++) {
    for (j=0; j<10; j++) {
      E = velocityVerlet_tip4p (sg, kc, klj, Fc, Flj, fc, flj, delta, t);
      time += delta;
    }
    E_rot = E_trans = 0.0;
    for (j=0; j<n; j++) {
      w = t->mol[j]->w;
      I = t->mol[j]->I;
      v = t->mol[j]->v;
      E_rot += 0.5 * (w[0] * w[0] * I[0] + w[1] * w[1] * I[1] 
            + w[2] * w[2] * I[2]);
      E_trans += 0.5 * t->mol[j]->m * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    }
    //printf ("t = %le  T = %lf  E_trans = %le  E_rot = %le  E_kin = %le\n", time, 2.0 * E / ((6 * n - 3) * K_B), E_trans, E_rot, E);
  }
  
  printf ("-------------------------------------------------------\n");
  printf ("Cleaning up\n");
  
  del_tip4p (t);
  del_spatialgeometry (sg);
  del_spatialcluster (sroot);
  del_block (broot);
  del_clusterbasis (cbc);
  del_clusterbasis (cblj);
  del_clusterbasis (rbc);
  del_clusterbasis (rblj);
  delParameters_coulomb (c);
  delParameters_lj (lj);
  del_kernelmatrix (kc);
  del_kernelmatrix (klj);
  del_h2matrix (Vc);
  del_h2matrix (Vlj);
  del_h2matrix (Fc);
  del_h2matrix (Flj);
  del_avector (flj);
  del_avector (fc);
  
  return EXIT_SUCCESS;
}