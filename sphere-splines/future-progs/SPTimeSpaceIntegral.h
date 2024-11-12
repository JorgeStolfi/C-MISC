/* SPTimeSpaceIntegral.h -- integration of functions over sphere × time */
/* Last edited on 2005-08-14 17:12:17 by stolfi */

#ifndef SPTimeSpaceIntegral_H
#define SPTimeSpaceIntegral_H

#include <SPTriang.h>
#include <vec.h>
#include <js.h>

/* CUSTOM INTEGRALS */

/* These procedures accumulate their result on the argument {*sum},
  and in the correction term {*corr}, using Kahan's summation
  formula. The client must initialize both {sum} and {corr} with 0. */

void SPTimeSpaceIntegral_BySamples
  ( double func(S2Point *p, double t),
    S2Point_vec_t sp,
    double_vec_t wp,
    double_vec_t st,
    double_vec_t wt,
    double *sum,
    double *corr
  );
  /* Computes {SUM { wp[k]*wt[j]*func(sp[k],st[j]) : k = 0..M-1, j=0..N-1 }},
    where {sp} is a client-specified set of sample points, {wp} their
    corresponding weights (measures of the domain elements), {st} is a
    set of sample times, and {wt} are the corresponding weights
    (interval widths). */
 
/* INTEGRAL ON THE SPHERE TIMES AN INTERVAL */

double SPTimeSpaceIntegral_OnSphereInterval
  ( double func(S2Point *p, double t),
    Triangulation *tri,
    int smpOrderTime,
    double tMin,
    double tMax
  );
  /* Computes the integral of {func(p)} for {p} ranging over the whole
    sphere and {t} ranging over the interval {[tMin _ tMax]}.
    
    The spherical integral uses the sample points and weights stored
    in {tri}. If {tri} is NULL, the procedure uses the default
    triangulation {SPIntegral_GetDefaultTriangulation()}. 
    
    The time integral uses Gaussian quadrature of order {smpOrderTime}. */

#endif
