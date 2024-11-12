/* SPTimeOnlyIntegral.h -- integration of functions over sphere × time */
/* Last edited on 2005-08-14 17:18:51 by stolfi */

#ifndef SPTimeOnlyIntegral_H
#define SPTimeOnlyIntegral_H

#include <vec.h>
#include <js.h>

/* SIMPLE INTEGRALS OF REAL FUNCTIONS OF A REAL VARIABLE */

/* These procedures accumulate their result on the argument {*sum},
  and in the correction term {*corr}, using Kahan's summation
  formula. The client must initialize both {sum} and {corr} with 0. */

void SPTimeOnlyIntegral_BySamples
  ( double func(double t),
    double_vec_t st,
    double_vec_t wt,
    double *sum,
    double *corr
  );
  /* Computes {SUM { wt[j]*func(st[j]) : j=0..N-1 }}, where {st} is a
    set of sample times, and {wt} are the corresponding weights
    (interval widths). */
 
/* INTEGRAL ON AN INTERVAL */

double SPTimeOnlyIntegral_OnInterval
  ( double func(double t),
    int smpOrder,
    double tMin,
    double tMax
  );
  /* Computes the integral of {func(t)} for {t} ranging over the
    interval {[tMin _ tMax]}. Uses Gaussian quadrature of order
    {smpOrder}. */

/* GENERATING SAMPLE POINTS IN AN INTERVAL */

/* For the tools in this section, the parameter {*ns} is the number of
  samples already stored in {st}. The new samples are stored starting
  at {st.e[*ns]}, and the value of {*ns} is incremented accordingly.
  The integration weight corresponding to {st.e[i]} is stored in
  {wt.e[i]}. The arrays {st} and {wt} are allocated and/or expanded
  as needed. */

void SPTimeOnlyIntegral_GaussSampleInterval
  ( double tMin, double tMax, 
    int smpOrder, 
    double_vec_t *st, 
    double_vec_t *wt, 
    int *ns
  );
  /* Stores in {st[*ns..*ns+smpOrder-1]} and
    {wt[*ns..*ns+smpOrder-1]} the knots and weights for Gaussian
    integration of order {smpOrder} over the time interval {[tMin _
    tMax]}. */
  
#endif
