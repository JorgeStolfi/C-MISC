/* SOIntegral.h -- integration of functions on the root cell */
/* Last edited on 2007-01-04 00:20:27 by stolfi */

#ifndef SOIntegral_H
#define SOIntegral_H

#include <SOGrid.h>

#include <dg_tree.h>
#include <dg_grid.h>

#include <vec.h>

/* INTEGRANDS */

typedef void SOIntegral_Func(double *x, double *fx);
  /* A generic function to be integrated: given {x[0..m-1]}, computes 
     {fx[0..n-1] = f(x)}.  We assume that the function itself knows 
     the dimension {m} of {x} and the dimension {n} of {fx}. */

/* SIMPLE INTEGRALS */

/* These procedures accumulate their result on the argument {sum[0..n-1]},
  and in the correction term {corr[0..n-1]}, using Kahan's summation
  formula. The client must initialize both {sum} and {corr} with 0. */

extern int SOIntegral_GaussOrder;
  /* Number of knots to be used by {SOIntegral_Gauss} below,
    along each coordinate axis. */

void SOIntegral_Gauss
  ( SOIntegral_Func f, 
    dg_dim_t d,     /* Dimension {m} of domain. */
    dg_dim_t fd,    /* Dimension {n} of the function's range. */
    double *sum, 
    double *corr
  );
  /* Computes the integral of {f(x)} over the {d}-dimensional box with 
    corners {(0,..)} and {(1,..)}, using Gaussian quadrature on a 
    sample grid with {SOIntegral_GaussOrder} knots along each axis.
    Accumulates the integral in {sum[0..fd-1]} and {corr[0..fd-1]}. */

void SOIntegral_Gauss_Limits
  ( SOIntegral_Func f, 
    dg_dim_t d, 
    dg_dim_t fd, 
    double *sum, 
    double *corr,
    double *min,    /* Lower coordinates of integration limits. */
    double *max     /* Higher coordinates of integration limits. */
  );

void SOIntegral_OnRootCell
  ( SOIntegral_Func f, 
    dg_dim_t d,     /* Dimension {m} of domain. */
    dg_dim_t fd,    /* Dimension {n} of the function's range. */
    SOGrid_Tree *tree,     /* Finite cell grid to use. */
    double *sum, 
    double *corr
  );
  /* Computes the integral of {f(x)} over the {d}-dimensional box with
    corners {(0,..)} and {(1,..)}, using {SOIntegral_Gauss} in each cell
    of the given {tree}. Accumulates the integral in {sum[0..fd-1]} and
    {corr[0..fd-1]}. */

#endif
