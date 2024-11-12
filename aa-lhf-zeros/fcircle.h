/* Last edited on 2010-08-15 00:45:41 by stolfi */
#include <aa.h>

#ifndef fcircle_H
#define fcircle_H

#ifndef real
#define	real	double
#endif

#define fcircle_N 2 /* number of variables */
#define fcircle_M 1 /* number of equations */
#define fcircle_P 2 /* Number of noises that may be eliminated */

void fcircle_func_args(real xmin, real xmax, real ymin, real ymax, AAform X[], int S[]);
  /* {X[0..N-1]}, {S[0..P-1]}. */

void fcircle_func_aa(AAform X[], AAform Y[]);
  /* {X[0..N-1]}, {Y[0..M-1]}. */

void fcircle_func_show(void);

#endif
