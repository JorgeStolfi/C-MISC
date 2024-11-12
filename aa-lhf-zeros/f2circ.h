/* Last edited on 2010-08-15 00:45:44 by stolfi */
#include <aa.h>

#ifndef f2circ_H
#define f2circ_H

#ifndef real
#define	real	double
#endif

#define f2circ_N 2 /* number of variables */
#define f2circ_M 2 /* number of equations */
#define f2circ_P 2 /* number of noises that may be eliminated */

void f2circ_func_args(real xmin, real xmax, real ymin, real ymax, AAform X[], int S[]);
  /* {X[0..N-1]}, {S[0..P-1]}. */

void f2circ_func_aa(AAform X[], AAform Y[]);
  /* {X[0..N-1]}, {Y[0..M-1]}. */

void f2circ_func_show(void);

#endif
