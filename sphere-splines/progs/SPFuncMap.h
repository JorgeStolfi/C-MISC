/* SPFuncMap.h -- pointwise operator for spherical functions. */
/* Last edited on 2005-06-05 14:27:05 by stolfi */

#ifndef SPFuncMap_H
#define SPFuncMap_H

#include <SPBasic.h>
#include <bool.h>

/* 
  An {SPFuncMap} is an object that describes a pointwise 
  operator {F(u,p)} that modifies the value {u} of a spherical
  function, possibly depending on the point {p}.
  
  This object is used mainly as the right-hand-side of the general
  differential equation {D(f)(p) == F(f(p),p)} where {D} is a
  linear spherical differential operator. It is also compatible with
  the {FuncMap} argument of {SPFunction_Dot} and its relatives.
*/

typedef double FuncMapMap(double u, S2Point *p); 
  /* A procedure that implements a pointwise spherical operator. */

typedef struct SPFuncMap   /* A pointwise spherical operator. */
  { FuncMapMap *map;   /* The map proper. */
    bool_t zeroPres;   /* TRUE implies zero-preserving, i.e. {F(0,p) = 0}. */
    char *descr;       /* The formula for {F(u, (x,y,z))}. */
  } SPFuncMap;

typedef SPFuncMap FuncMap; /* A shorter name. */

/* The identity operator: */
#define NoFMap (FuncMap){ NULL, TRUE, "u" }  

SPFuncMap SPFuncMap_FromName(char *name);
  /* Returns a {FuncMap} and its data given its name.
    Bombs out if the name is unrecognized. */

#endif

