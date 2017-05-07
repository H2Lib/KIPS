
#include <stdio.h>
#include "interpolation.h"
#include "basic.h"

#ifdef USE_FLOAT
static const real tolerance = 5.0e-5;
#else
static const real tolerance = 1.0e-12;
#endif

int
main(int argc, char **argv)
{
  pinterpolation in1, in2;
  uint m;
  real a, b;
  real error;
  uint problems = 0;
  uint i;

  init_kips(&argc, &argv);

  m = 19;
  
  (void) printf("Creating Chebyshev points of degree %u\n",
		m-1);
  in1 = new_interpolation(m);
  chebyshev_interpolation(-1.0, 1.0, in1);

  error = 0.0;
  for(i=0; i<m; i++)
    error = REAL_MAX(error,
		     REAL_ABS(in1->xi[i] - REAL_COS(M_PI * (i+0.5) / m)));
  (void) printf("  Error %.3e",
		error);
  if(error <= tolerance)
    (void) printf("  --  Okay\n");
  else {
    (void) printf("  --  NOT Okay\a\n");
    problems++;
  }

  a = -(REAL_RAND() + 1.5);
  b = REAL_RAND() + 1.5;

  (void) printf("Transforming to [%.3f,%.3f]\n",
		a, b);
  in2 = build_transformed_interpolation(in1, a, b);
  chebyshev_interpolation(a, b, in1);

  error = 0.0;
  for(i=0; i<m; i++)
    error = REAL_MAX(error,
		     REAL_ABS(in1->xi[i] - in2->xi[i]));
  (void) printf("  Error %.3e",
		error);
  if(error <= tolerance)
    (void) printf("  --  Okay\n");
  else {
    (void) printf("  --  NOT Okay\a\n");
    problems++;
  }


  (void) printf("----------------------------------------\n"
		"%u errors found\n",
		problems);

  uninit_kips();

  return 0;
}
