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
  preal bmin, bmax, E_ref, *xn, *xc, r, com;
  uint dim, i, m, maxdepth, n;
  real eta, maxdiam, R, E, sig, eps, alpha, emax;
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
  char *file;
  pavector fc, flj;
  pquaternion q = new_quaternion ();
  
  init_kips (&argc, &argv);
  
  dim = 3;
  maxdepth = 0;
  maxdiam = 1.6;
  alpha = 0.0;
  m = 3;
  eta = 1.0;
  sig = SIG_TIP4P;
  eps = EPS_TIP4P;
  s = &split_cutoff;
  sd = &derivative_cutoff;
  R = 20.0;
  file = (char *) allocmem (sizeof (char) * 100);
  bmin = allocreal (dim);
  bmax = allocreal (dim);
  E_ref = allocreal (20);
  E_ref[0] = -26087.57;
  E_ref[1] = -69993.87;
  E_ref[2] = -116590.42;
  E_ref[3] = -152109.00;
  E_ref[4] = -197780.53;
  E_ref[5] = -243572.40;
  E_ref[6] = -305518.32;
  E_ref[7] = -344435.21;
  E_ref[8] = -391022.68;
  E_ref[9] = -431489.07;
  E_ref[10] = -492908.24;
  E_ref[11] = -532970.54;
  E_ref[12] = -582991.10;
  E_ref[13] = -628372.25;
  E_ref[14] = -681192.59;
  E_ref[15] = -723807.75;
  E_ref[16] = -773231.87;
  E_ref[17] = -821036.92;
  E_ref[18] = -872988.75;
  E_ref[19] = -916707.43;
  
  printf ("Setting up bounding box and cluster and block trees\n");
  for (i=0; i<dim; i++) {
    bmin[i] = -20.0;
    bmax[i] = 20.0;
  }
  
  sg = new_spatialgeometry (dim, bmin, bmax);
  sroot = init_spatialgeometry (maxdepth, maxdiam, sg);
  c = setparametersGradient_coulomb (s, sd, alpha, R, sg);
  lj = setparameters_lj (R, sig, eps, sg);
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

  printf ("Setting up H2-matrix structures");

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
  /*
  printf ("Calculating potential energies of water clusters\n");
  for (i=2; i<=21; i++) {
    printf ("%d molecules\n", i); 
    sprintf (file, "/home/jol/Dokumente/KIPS/trunk/Tests/TIP4P-%d.xyz", i);
    t = inputfile_tip4p (file);
    E = N_A * energy_tip4p (sg, kc, klj, Vc, Vlj, t);
    printf ("E = %16.14e J/mol\n", E);
    printf ("Relative error: %e\n", REAL_ABS ((E-E_ref[i-2])/E_ref[i-2]));
    printf ("\n");
    
    del_tip4p (t);
  }
  
  printf ("-------------------------------------------------------\n");
  printf ("Calculating forces of water clusters\n");
  
  for (i=2; i<=21; i++) {
    printf ("%d molecules\n", i); 
    sprintf (file, "/home/jol/Dokumente/KIPS/trunk/Tests/TIP4P-%d.xyz", i);
    t = inputfile_tip4p (file);
    n = t->n;
    fc = new_avector (3 * n * 3);
    flj = new_avector (n * 3);
    
    update_kernelmatrix (3 * n, t->xc, t->xcmol, kc);
    fullForce_coulomb (t->q, kc, sg, Fc, fc);
    
    update_kernelmatrix (n, t->xn, t->xnmol, klj);
    fullForce_lj (klj, sg, Flj, flj);
    
    emax = calcForce_tip4p (fc, flj, t);
    printf ("Maximal relative error of molecular net forces: %e\n", emax);
    emax = calcTorque_tip4p (fc, flj, t);
    printf ("Maximal relative error of molecular net torques: %e\n", emax);
    
    del_avector (fc);
    del_avector (flj);
    del_tip4p (t);
  }
  
  printf ("-------------------------------------------------------\n");
  printf ("Optimization via gradient method for two molecules\n");
  
  sprintf (file, "/home/jol/Dokumente/KIPS/trunk/Tests/TIP4P-2b.xyz");
  t = inputfile_tip4p (file);
  n = t->n;
  fc = new_avector (3 * n * 3);
  flj = new_avector (n * 3);

  update_kernelmatrix (3 * n, t->xc, t->xcmol, kc);
  fullForce_coulomb (t->q, kc, sg, Fc, fc);

  update_kernelmatrix (n, t->xn, t->xnmol, klj);
  fullForce_lj (klj, sg, Flj, flj);
  
  E = N_A * energy_tip4p (sg, kc, klj, Vc, Vlj, t);
  printf ("E = %16.14e J/mol\n", E);
  printf ("Relative error: %e\n", REAL_ABS ((E-E_ref[0])/E_ref[0]));

  emax = calcForce_tip4p (fc, flj, t);
  printf ("Maximal relative error of molecular net forces: %e\n", emax);
  emax = calcTorque_tip4p (fc, flj, t);
  printf ("Maximal relative error of molecular net torques: %e\n", emax);
  
  printf ("\n");
  printf ("Applying 100 steps of gradient descent method\n");
  for (i=0; i<100; i++) {
    E = N_A * gradientDescent_tip4p (sg, kc, klj, Vc, Fc, Vlj, Flj, 1.0e20, 0.5, 100, t);
  }
  

  printf ("E_new = %16.14e J/mol\n", E);
  printf ("Relative error: %e\n", REAL_ABS ((E-E_ref[0])/E_ref[0]));
  printf ("\n");
  
  update_kernelmatrix (3 * n, t->xc, t->xcmol, kc);
  fullForce_coulomb (t->q, kc, sg, Fc, fc);

  update_kernelmatrix (n, t->xn, t->xnmol, klj);
  fullForce_lj (klj, sg, Flj, flj);
  
  emax = calcForce_tip4p (fc, flj, t);
  printf ("Maximal relative error of molecular net forces: %e\n", emax);
  emax = calcTorque_tip4p (fc, flj, t);
  printf ("Maximal relative error of molecular net torques: %e\n", emax);
  
  del_avector (fc);
  del_avector (flj);
  del_tip4p (t);
  */
  printf ("-------------------------------------------------------\n");
  printf ("Optimizing orientation\n");
  
  sprintf (file, "/home/jol/Dokumente/KIPS/trunk/Tests/TIP4P-2.xyz");
  t = inputfile_tip4p (file);
  n = t->n;
  fc = new_avector (3 * n * 3);
  flj = new_avector (n * 3);
  
  q[0] = 0.8;
  q[1] = 0.6;
  q[2] = 0.0;
  q[3] = 0.0;
  inertia_tip4p (t);
  refOrientation_molecule (t->mol[0]);
  adjust_molecule (q, t->mol[0]);
  xn = t->xn;
  xc = t->xc;
  com = t->mol[0]->com->v;
  r = t->mol[0]->r[0]->v;
  xn[0][0] = com[0] + r[0];
  xn[0][1] = com[1] + r[1];
  xn[0][2] = com[2] + r[2];
  r = t->mol[0]->r[1]->v;
  xc[1][0] = com[0] + r[0];
  xc[1][1] = com[1] + r[1];
  xc[1][2] = com[2] + r[2];
  r = t->mol[0]->r[2]->v;
  xc[2][0] = com[0] + r[0];
  xc[2][1] = com[1] + r[1];
  xc[2][2] = com[2] + r[2];
  calc_virtual_charge (xc[1], xc[2], xn[0], xc[0]);

  update_kernelmatrix (3 * n, t->xc, t->xcmol, kc);
  fullForce_coulomb (t->q, kc, sg, Fc, fc);

  update_kernelmatrix (n, t->xn, t->xnmol, klj);
  fullForce_lj (klj, sg, Flj, flj);
  
  E = N_A * energy_tip4p (sg, kc, klj, Vc, Vlj, t);
  printf ("E = %16.14e J/mol\n", E);
  printf ("Relative error: %e\n", REAL_ABS ((E-E_ref[0])/E_ref[0]));

  emax = calcForce_tip4p (fc, flj, t);
  printf ("Maximal relative error of molecular net forces: %e\n", emax);
  emax = calcTorque_tip4p (fc, flj, t);
  printf ("Maximal relative error of molecular net torques: %e\n", emax);
  
  printf ("\n");
  printf ("Applying 100 steps of gradient descent method\n");
  for (i=0; i<100; i++) {
    //E = N_A * gradientDescent_tip4p (sg, kc, klj, Vc, Fc, Vlj, Flj, 1.0e20, 0.5, 100, t);
    E = N_A * gradientDescentRotation_tip4p (sg, kc, klj, Vc, Fc, Vlj, Flj, 1.0e30, 0.5, 100, t);
  }
  
  
  

  printf ("E_new = %16.14e J/mol\n", E);
  printf ("Relative error: %e\n", REAL_ABS ((E-E_ref[0])/E_ref[0]));
  printf ("\n");
  
  update_kernelmatrix (3 * n, t->xc, t->xcmol, kc);
  fullForce_coulomb (t->q, kc, sg, Fc, fc);

  update_kernelmatrix (n, t->xn, t->xnmol, klj);
  fullForce_lj (klj, sg, Flj, flj);
  
  emax = calcForce_tip4p (fc, flj, t);
  printf ("Maximal relative error of molecular net forces: %e\n", emax);
  emax = calcTorque_tip4p (fc, flj, t);
  printf ("Maximal relative error of molecular net torques: %e\n", emax);
  
  del_avector (fc);
  del_avector (flj);
  del_tip4p (t);
  
  printf ("-------------------------------------------------------\n");
  printf ("Cleaning up\n");
  del_spatialgeometry (sg);
  del_spatialcluster (sroot);
  del_block (broot);
  del_clusterbasis (cbc);
  del_clusterbasis (cblj);
  del_clusterbasis (rbc);
  del_clusterbasis (rblj);
  freemem (c);
  freemem (lj);
  del_kernelmatrix (kc);
  del_kernelmatrix (klj);
  del_h2matrix (Vc);
  del_h2matrix (Vlj);
  del_h2matrix (Fc);
  del_h2matrix (Flj);
  freemem (file);
  freemem (E_ref);
  del_quaternion (q);
  
  return EXIT_SUCCESS;
}