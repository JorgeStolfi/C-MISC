/* SPTimeOnlyFuncMap.h -- pointwise operator for sphere-time functions. */
/* Last edited on 2005-08-14 17:30:05 by stolfi */

#ifndef SPTimeOnlyFuncMap_H
#define SPTimeOnlyFuncMap_H

#include <SPBasic.h>
#include <bool.h>

/* 
  An {SPTimeOnlyFuncMap} is an object that describes a pointwise 
  operator {F(u,t)} that modifies the value {u} of a time function
  function, possibly depending on the time {t}.
  
  This object is used mainly as the right-hand-side of the general
  differential equation {D(f)(t) == F(f(t),t)} where {D} is a
  linear sphere-time differential operator. It is also compatible with
  the {TimeOnlyFuncMap} argument of {SPTimeOnlyFunction_Dot}.
*/

typedef double TimeOnlyFuncMapMap(double u, double t); 
  /* A procedure that implements a pointwise sphere-time operator. */

typedef struct SPTimeOnlyFuncMap   /* A pointwise sphere-time operator. */
  { TimeOnlyFuncMapMap *map;   /* The map proper. */
    bool_t zeroPres;    /* TRUE implies zero-preserving, i.e. {F(0,t) = 0}. */
    char *descr;        /* The formula for {F(u,t)}. */
  } SPTimeOnlyFuncMap;
  /* A pointwise operator that sphere-time suitable  */

typedef SPTimeOnlyFuncMap TimeOnlyFuncMap; /* A shorter name. */

/* The identity operator: */
#define NoTOFMap (TimeOnlyFuncMap){ NULL, TRUE, "u" }  

SPTimeOnlyFuncMap SPTimeOnlyFuncMap_FromName(char *name);
  /* Returns a {SPTimeOnlyFuncMap} and its data given its name.
    Bombs out if the name is unrecognized. */

#endif

