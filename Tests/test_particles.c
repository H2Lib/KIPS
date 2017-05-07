
#include <stdio.h>

#include "particles.h"
#include "cluster.h"
#include "block.h"

int
main(int argc, char **argv)
{
  pparticles p;
  pclustergeometry cg;
  uint *idx;
  pcluster root;
  pblock broot;
  pavector m, phi;
  uint n;
  field alpha;
  uint i;

  init_kips(&argc, &argv);

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

  uninit_kips();

  return 0;
}
