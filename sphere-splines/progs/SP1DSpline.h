/* SP1DSpline.h -- unidimensional splines. */
/* Last edited on 2003-12-21 18:57:55 by stolfi */

#ifndef SP1DSpline_H
#define SP1DSpline_H

/*
  This module defines finite-element bases for univariate splines 
  over the integer grid.  
  
  The basis elements have a specified continuity class {c} and degree
  {g = 2·c + 1}. There are {c + 1} basis elements associated with each
  knot {j}; they are all translates of the {c+1} elements
  {bas[0..c](t)} associated with knot {0}. Their support consists of
  the two unit intervals adjacent to that knot, i.e. the interval
  {[j-1 _ j+1]}.
  
  The elements are /Hermite-like/ in the sense that the value and all
  derivatives to order {c} of element {bas[i]} are equal to zero,
  except for the derivative of order {i}. Elements with even {i} are
  symmetric around the central knot, elements with odd {i} are
  anti-symmetric. */
  
double SP1DSpline_Eval(int c, int i, double t); 
  /* Computes the Hermite-like spline element {bas[i](t)} of
    continuity {c} and index {i} (which must lie in {0,..c}) for the
    argument {t}.  */

double SP1DSpline_Eval_0(int c, int i); 
  /* Same as {SP1DSpline_Eval(c,i,0)}, but probably faster. */

void SP1DSpline_Coeffs(int c, int i, int j, double *a);
  /* Returns in {a[0..d]} the coefficients of the Hermite-like spline
    element {bas[i](t)] of continuity {c} and index {i}, restricted to
    the interval {[j..j+1]}, in terms of the variable {z = t-j}. The
    coefficients are such that {bas[i](j+z) = SUM { a[i]·z^i : i =
    0..d}} where {z = t} ranges in {[0 _ 1]}, and {d = 2·c + 1}. */

#endif
