/* SOFuncMap.h -- right-hand-side functions for diff eqns. */
/* Last edited on 2007-01-04 00:22:09 by stolfi */

#ifndef SOFuncMap_H
#define SOFuncMap_H

#include <SOFunction.h>
#include <SOBasic.h>

#include <dg_grid.h>

typedef FuncMap SOFuncMap;
  /* A function suitable for the right-hand-side of the general
    differential equation {D(f)(p) == F(f(p),p)} where {D} is a
    differential operator.  Also compatible with the {FuncMap}
    argument of {SOFunction_Dot}. */

typedef struct SOFuncMap_Data 
  { FuncMap FMap;   /* The function map. */
    char *name;     /* The identifying name for {FMap}. */
    char *descr;    /* The formula for {FMap(u, (x,y,z))}. */
    double coeff;   /* Coefficient for the screened Poisson equation. */
    char *sol;      /* screened Poisson solution for {coeff} and {FMap} as RHS. */
  } SOFuncMap_Data;
  /* A function map {FMap} and some data about it.
    The {sol} string, if not NULL, is the name of an 
    {SOProcFunction} {f} which is the exact solution 
    of the screened Poisson equation 
    {Lap(f)(p) + c*f(p) == FMap(f(p),p)},
    when {c == coeff}. */

SOFuncMap_Data SOFuncMap_FromName(char *name, dg_dim_t u, dg_dim_t p, dg_dim_t v);
  /* Returns a {FuncMap} and its data given its name.
    Bombs out if the name is unrecognized. */

#endif

