/* Last edited on 2010-08-15 00:45:37 by stolfi */
#include <aa.h>
#ifndef fcirhex_H
#define fcirhex_H

#ifndef real
#define	real	double
#endif

#define fcirhex_N 2 /* number of variables */
#define fcirhex_M 1 /* number of equations */
#define fcirhex_P 3 /* Number of noises that may be eliminated */

void fcirhex_func_args(real xmin, real xmax, real ymin, real ymax, AAform X[], int S[]);
  /* {X[0..N-1]}, {S[0..P-1]}. */

void fcirhex_func_aa(AAform X[], AAform Y[]);
  /* {X[0..N-1]}, {Y[0..M-1]}. */

void fcirhex_func_show(void);

#endif
