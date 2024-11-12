/* Time elements with degree 3 continuity 1. */
/* Last edited on 2005-08-27 21:18:29 by stolfi */

#ifndef SPTimeBasisG3C1_H
#define SPTimeBasisG3C1_H

#include <bool.h>

/* In the following procedures, {TD(f,diff)} means the {diff}th time derivative
   of function {f}; {tg} is the maximum degree of a time element (3); and 
   {pg} is the maximum degree of a gauge element (1). */

double SPTimeBasisG3C1_BasisEval(int u, double t, double tj, double tStep, int diff);
  /* Evaluates {TD(tau[j,u],diff)(t)}, i.e. {TD(tau[0,u],diff)(t - tj)}. */

void SPTimeBasisG3C1_BasisCoeffs(int k, int u, double tStep, int diff, double *tc);
  /* Returns in {tc[0..tg-diff]} the coefficients of {TD(tau[j-k,u],diff)(t)} 
    for {t} in the interval {[tj-tStep _ tj]}, that is,
    of {TD(tau[0,u],diff)(t - tj + k*tStep)} 
    
    The coefficients are such that 
      {TD(tau[j-k,u],diff)(t) = SUM { tc[i]·z^i : i = 0..tg-diff}} 
    where {z = (t-tj)/tStep+1} ranges in {[0 _ 1]}. Note that the array {tc} 
    must have {tg+1} elements, even though only {tg-diff+1} are needed 
    for the result. */

double SPTimeBasisG3C1_GaugeEval(int v, double t, double tj, double tStep);
  /* Evaluates {pi[0,v](t)}. */

void SPTimeBasisG3C1_GaugeCoeffs(int v, double *pc);
  /* Returns in {pc[0..pg]} the coefficients of the time gauge element
    {pi[j,u]} in the time interval {[j-1 _ j]}.  The coefficients 
    are are such that {pi[j,u] = SUM { pc[i]·z^i : i = 0..1}}
    where {z = t-j+1} ranges in {[0 _ 1]}. */
    
double SPTimeBasisG3C1_BasisGaugeDot
  ( int k, 
    int u, 
    int v, 
    double tStep, 
    int diff, 
    bool_t verbose
  );
  /* These procedures compute the time-wise scalar products {
    <TD(tau[0,u],diff) | pi[k,v]> }, where {tau[0,u]} is the standard
    time basis element, {pi} is the gauge basis element,
    {TD(func,diff)} is the time derivative of {func} of order {diff}.
    Note that the result is independent of the epoch {j}. */

#endif

