
#include <stdio.h>

#include "particles.h"
#include "parameters.h"
#include "cluster.h"
#include "block.h"
#include "h2matrix.h"

int
main(int argc, char **argv)
{
  pparticles p;
  pclustergeometry cg;
  uint *idx;
  pcluster root;
  pblock broot;
  pclusterbasis cb;
  ph2matrix G;
  pavector m, phi;
  pstopwatch sw;
  real t_run;
  size_t sz;
  uint n;
  field alpha;
  uint i;

  init_kips(&argc, &argv);

  sw = new_stopwatch();

  n = askforint("Number of particles?", "kips_particles", 40000);

  (void) printf("Creating %u random particles\n",
		n);
  p = new_particles(n);
  random_particles(p);

  (void) printf("Setting up cluster tree\n");
  cg = buildgeometry_particles(p);
  idx = (uint *) allocmem(sizeof(uint) * n);
  for(i=0; i<n; i++)
    idx[i] = i;
  root = buildgeometric_clustergeometry(cg, 100, 64, n, idx);
  (void) printf("  %u clusters\n"
		"  %u levels\n",
		root->desc,
		getdepth_cluster(root) + 1);

  (void) printf("Setting up block tree\n");
  broot = buildh2std_block(root, root, 2.0);
  (void) printf("  %u blocks\n"
		"  %u levels\n",
		broot->desc,
		getdepth_block(broot) + 1);

  (void) printf("Setting up Lagrange cluster basis\n");
  start_stopwatch(sw);
  cb = build_fromcluster_clusterbasis(root);
  t_run = stop_stopwatch(sw);
  sz = getsize_clusterbasis(cb);
  (void) printf("  %.1f seconds\n"
		"  %.2f MB (%.1f KB/DoF)\n",
		t_run,
		sz / 1048576.0, sz / 1024.0 / n);

  (void) printf("Setting up H2-matrix\n");
  start_stopwatch(sw);
  G = build_fromblock_h2matrix(broot, cb, cb);
  t_run = stop_stopwatch(sw);
  sz = getsize_h2matrix(G);
  (void) printf("  %.1f seconds\n"
		"  %.2f MB (%.1f KB/DoF)\n",
		t_run,
		sz / 1048576.0, sz / 1024.0 / n);

  uninit_kips();

  return 0;
}
