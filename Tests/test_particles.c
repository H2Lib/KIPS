
#include <stdio.h>

#include "particles.h"
#include "parameters.h"
#include "cluster.h"
#include "block.h"
#include "h2matrix.h"

#ifdef USE_FLOAT
static const real tolerance = 5.0e-5;
#else
static const real tolerance = 1.0e-12;
#endif

int
main(int argc, char **argv)
{
  pparticles p;
  pclustergeometry cg;
  uint     *idx;
  pcluster  root;
  pblock    broot;
  pclusterbasis cb;
  ph2matrix G;
  pavector  src, trg;
  pamatrix  V, V1, Vson, E;
  pstopwatch sw;
  real      norm, error;
  real      t_run;
  uint      problems = 0;
  size_t    sz;
  uint      n, m, res;
  field     alpha;
  uint      i, off;

  init_kips(&argc, &argv);

  sw = new_stopwatch();

  n = askforint("Number of particles?", "kips_particles", 40000);

  m = askforint("Interpolation order?", "kips_intorder", 4);

  (void) printf("Creating %u random particles, interpolation order %u\n",
		n, m);
  p = new_particles(n, m);
  random_particles(p);

  res = 2 * m * m * m;
  (void) printf("Setting up cluster tree, resolution %u\n", res);
  cg = buildgeometry_particles(p);
  idx = (uint *) allocmem(sizeof(uint) * n);
  for (i = 0; i < n; i++)
    idx[i] = i;
  root = buildgeometric_clustergeometry(cg, 100, 2 * m * m * m, n, idx);
  (void) printf("  %u clusters\n"
		"  %u levels\n", root->desc, getdepth_cluster(root) + 1);

  (void) printf("Creating V matrix for the root\n");
  V = new_amatrix(root->size, m * m * m);
  buildV_particles(p, root, V);
  norm = normfrob_amatrix(V);

  if (root->sons > 0) {
    (void) printf("Subtracting V E for the sons\n");
    off = 0;
    for (i = 0; i < root->sons; i++) {
      V1 = new_sub_amatrix(V, root->son[i]->size, off, m * m * m, 0);
      Vson = new_amatrix(root->son[i]->size, m * m * m);
      buildV_particles(p, root->son[i], Vson);

      E = new_amatrix(m * m * m, m * m * m);
      buildE_particles(p, root->son[i], root, E);
      addmul_amatrix(-1.0, false, Vson, false, E, V1);

      del_amatrix(E);
      del_amatrix(Vson);
      del_amatrix(V1);

      off += root->son[i]->size;
    }
    assert(off == root->size);

    error = normfrob_amatrix(V);
    (void) printf("  Error %.3e (%.3e)", error, error / norm);
    if (error <= tolerance * norm)
      (void) printf("  --  Okay\n");
    else {
      (void) printf("  --  NOT Okay\a\n");
      problems++;
    }
  }
  del_amatrix(V);

  (void) printf("Setting up block tree\n");
  broot = buildh2std_block(root, root, 2.0);
  (void) printf("  %u blocks\n"
		"  %u levels\n", broot->desc, getdepth_block(broot) + 1);

  (void) printf("Setting up Lagrange cluster basis\n");
  start_stopwatch(sw);
  cb = build_fromcluster_clusterbasis(root);
  fill_clusterbasis((buildV_t) buildV_particles,
		    (buildE_t) buildE_particles, p, cb);
  t_run = stop_stopwatch(sw);
  sz = getsize_clusterbasis(cb);
  (void) printf("  %.1f seconds\n"
		"  %.2f MB (%.1f KB/DoF)\n"
		"  k %u, ktree %u\n",
		t_run, sz / 1048576.0, sz / 1024.0 / n, cb->k, cb->ktree);

  (void) printf("Setting up H2-matrix\n");
  start_stopwatch(sw);
  G = build_fromblock_h2matrix(broot, cb, cb);
  fill_h2matrix((buildN_t) buildN_particles,
		(buildS_t) buildS_particles, p, G);
  t_run = stop_stopwatch(sw);
  sz = getsize_h2matrix(G);
  (void) printf("  %.1f seconds\n"
		"  %.2f MB (%.1f KB/DoF)\n",
		t_run, sz / 1048576.0, sz / 1024.0 / n);

  (void) printf("Testing approximation\n");
  src = new_avector(n);
  trg = new_avector(n);
  random_avector(src);
  alpha = FIELD_RAND();
  clear_avector(trg);

  (void) printf("  Direct evaluation\n");
  start_stopwatch(sw);
  addeval_direct_particles(alpha, p, src, trg);
  t_run = stop_stopwatch(sw);
  norm = norm2_avector(trg);
  (void) printf("  %.2f seconds\n" "  Norm %.3e\n", t_run, norm);

  (void) printf("  H^2-matrix evaluation\n");
  start_stopwatch(sw);
  addeval_h2matrix(-alpha, G, src, trg);
  t_run = stop_stopwatch(sw);
  error = norm2_avector(trg);
  (void) printf("  %.2f seconds\n"
		"  Error %.3e (%.3e)", t_run, error, error / norm);

  if (error <= 15.0 * pow(0.125, m) * norm)
    (void) printf("  --  Okay\n");
  else {
    (void) printf("  --  NOT Okay\a\n");
    problems++;
  }

  (void) printf("Testing on-the-fly H^2-matrix evaluation\n");
  random_avector(src);
  clear_avector(trg);
  start_stopwatch(sw);
  addeval_h2matrix(alpha, G, src, trg);
  t_run = stop_stopwatch(sw);
  norm = norm2_avector(trg);
  (void) printf("  %.2f seconds\n"
		"  Norm %.3e\n",
		t_run, norm);
  start_stopwatch(sw);
  addeval_otf_h2matrix(-alpha, broot, cb, cb,
		       (buildN_t) buildN_particles,
		       (buildS_t) buildS_particles, p, src, trg);
  t_run = stop_stopwatch(sw);
  error = norm2_avector(trg);
  (void) printf("  %.2f seconds\n"
		"  Error %.3e (%.3e)",
		t_run, error, error / norm);
  if (error <= tolerance * norm)
    (void) printf("  --  Okay\n");
  else {
    (void) printf("  --  NOT Okay\a\n");
    problems++;
  }

  (void) printf("Testing on-the-fly adjoint H^2-matrix evaluation\n");
  random_avector(src);
  clear_avector(trg);
  start_stopwatch(sw);
  addevaltrans_h2matrix(alpha, G, src, trg);
  t_run = stop_stopwatch(sw);
  norm = norm2_avector(trg);
  (void) printf("  %.2f seconds\n"
		"  Norm %.3e\n",
		t_run, norm);
  start_stopwatch(sw);
  addevaltrans_otf_h2matrix(-alpha, broot, cb, cb,
			    (buildN_t) buildN_particles,
			    (buildS_t) buildS_particles, p, src, trg);
  t_run = stop_stopwatch(sw);
  error = norm2_avector(trg);
  (void) printf("  %.2f seconds\n"
		"  Error %.3e (%.3e)",
		t_run, error, error / norm);
  if (error <= tolerance * norm)
    (void) printf("  --  Okay\n");
  else {
    (void) printf("  --  NOT Okay\a\n");
    problems++;
  }

  del_avector(trg);
  del_avector(src);
  del_h2matrix(G);
  del_clusterbasis(cb);

  (void) printf("----------------------------------------\n"
		"  %u matrices and\n"
		"  %u vectors still active\n"
		"  %u errors found\n",
		getactives_amatrix(), getactives_avector(), problems);

  uninit_kips();

  return 0;
}
