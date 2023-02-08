/* SOBasic.h -- basic data types for the Simoleo project. */
/* Last edited on 2008-07-14 22:22:23 by stolfi */

#ifndef SOBasic_H
#define SOBasic_H

#include <SOGrid.h>

#include <dg_grid.h>

#include <values.h>

typedef enum { NEG = -1, ZER = 0, POS = 1 } sign;

typedef char * string;

typedef double ScalarField (dg_dim_t m, double *x);
  /* A generic function from {R^m} to {R}: given {x[0..m-1]}, returns {f(x)}. */

typedef void VectorField (dg_dim_t m, dg_dim_t n, double *x, double *y);
  /* A generic function from {R^m} to {R^n}: given {x[0..m-1]},
    stores {f(x)} in {y[0..n-1]}. */

typedef enum{X, Y, Z, T} dimensions;

#define INFTY    INFINITY
    
#define TWOPI    (2*M_PI)
#define FOURPI   (4*M_PI)
#define PIBYFOUR (M_PI/4)

#define PHI      (1.6180339887498948482)
#define SQRT3    (1.73205080756887729352)
#define SQRTHALF (0.70710678118654752440084)

#define MAX_PDIM (dg_MAX_GRID_DIM)
  /* Maximum dimension {m} for the domain of an {SOFunction}. High dimensions 
    are dangerous because too many things have size {2^d}. */

#define MAX_FDIM (16)
  /* Maximum dimension {n} for the range space of an {SOFunction}. 
    High dimensions are not so dangerous here, except for splines. */

typedef struct Color { double c[3]; } Color;
  /* RGB intensities in {[0_1]}. */

/* Some useful macros: */

#define mumble(level, ...) \
  do { if (verbose > (level)) { fprintf(stderr, __VA_ARGS__); } } while(0)

#endif
